# M5 plan (draft, 2026-10-22)

M5 as listed: CFADRIZ, grid.reference = 'wcs', in-place append (§10.3), band-by-band history I/O, v0.3 → v0.4 converter.

This plan is based on the spec as it stands (v0.4 draft, uncommitted edits included) and on the code in src/lsf. Section A covers the spec review, section B the work items, section C the decisions needed before starting.

## A. Spec review: gaps and problems

Grouped by how much they matter. "Code" means the spec is right and the code is behind. "Spec" means the text itself needs to change.

### A1. Problems to fix in this wave

1. **§6 omits the normalization baseline (spec).** The writer stores `norm.estimator`, `norm.rgb_equal` and `norm.baseline` (the frozen loc/scale of the first batch's reference frame, lsf_settings.c). Appends depend on them. §6 lists only `norm.method`. Without `norm.baseline` an append cannot normalize the way the first batch did, so it is a load-bearing fact, not a detail. Add all three to the required keys, with `norm.baseline` conditional on norm.method ≠ none.
2. **LSF_REJ is in the spec but not in the code, and the code has a latent bug (code).** §5.4 describes the HDU and an option on lsf_init/lsf_append; neither exists. lsf_maps warns "not supported". On rewrite, lsf_fits.c drops an existing LSF_REJ but writes `rejmaps.enabled` back unchanged. A file with rejmaps = T would therefore fail check 7 on its next open, so the writer would corrupt its own output. No current writer produces LSF_REJ, so nobody hits this today. The per-pixel low/high counts are already computed in the merge loop (`rej[2]` from `apply_rejection_float`, plus the history test), so the feature is cheap. See C3.
3. **§10.3 in-place append has more problems than it solves (spec).**
   - It saves almost no I/O. An append changes MEAN/W/M2/W2/COUNT wherever the batch lands, which is nearly the whole grid. So in place we still rewrite essentially all of LSF_ACC and LSF_COUNT plus HDU 0. What it actually saves is disk space (no temporary copy) and the Windows rename step.
   - It is not crash-safe. LSFBUSY only detects an interruption, and recovery then needs a backup, which costs the disk space in-place was meant to save. The spec itself says a writer "SHOULD keep a backup".
   - Header growth. Each append adds HISTORY cards to HDU 0. Once the header crosses a 2880-byte block boundary, CFITSIO inserts a block and shifts every following byte of a multi-GB file. In-place would need reserved header space (blank cards written at init), and the spec does not say so.
   - It is incompatible with compression and with grid expansion. Both force a full rewrite anyway.

   Recommendation: drop §10.3 from M5. Either delete it or move it to §14 with the header-reservation requirement written down (C1).
4. **§10.11 converter claims "without loss" (spec).** Inverting output normalization (MEAN = p·(max−min)+min) in float32 is not bit-exact, because each pixel can be off by about 1 ulp. A v0.3 file that was output-normalized therefore converts with rounding error. The section should say "near-lossless (≤ 1 ulp) when LSFNORM = T". Separately: v0.3 files were only ever written by the Python POC (lsf_poc/lsf_fits.py, PROFILE_VERSION 0.3), and the C writer started at 0.4. So the converter is only worth building if there are v0.3 stacks you want to keep (C2).
5. **Facts versus policy in §6 (spec).** The agreed rule asks whether a reader needs the value to interpret the data. M5 adds new settings, and three of them fail that test:
   - `grid.reference` ('stack' / 'wcs') is the method of the next registration. It changes nothing about how the file is read, so it should be an lsf_append option, not a stored row (C4).
   - `grid.pixfrac` and `grid.kernel` are drizzle choices per append. Whatever a batch used should go into its ledger rows as a fact, not into a frozen setting.
   - `grid.oversample` stays in settings. It is the grid scale, which is a fact.
   - Existing rows (`align.*`, `reject.*`, `reject.hist_*`) act as recalled defaults, which §6 already says ("Settings are recalled on append"). They can stay. One sentence in §6 should separate "facts the data depends on" (sensor, cfa.mode, oversample, origin, weight constants, norm baseline, outnorm) from "defaults recalled on append".
6. **§10.6 memory refusal (spec).** "Refuses one whose expanded accumulator would not fit in the memory it may use" stops being true once band-by-band lands. The rule becomes: refuse only if one band of the expanded grid plus the batch buffers does not fit.
7. **The memory budget in lsf_process_batch ignores the history (code).** `max_rows` is computed from the full `get_max_memory_in_MB()` for frame buffers alone, but the whole history (18 B/px/chan) is already allocated on top of that. On a large grid the process can therefore exceed Siril's memory limit. This is fixed naturally by the band-by-band budget (B3).

### A2. Gaps that M5 features expose

8. **CFADRIZ per-pixel weights are not defined anywhere (spec, §10.4 / §10.8).** In Siril's drizzle stack the sample weight is frame weight × drizzle pixel weight (`dstack`, median_and_mean.c), and a sample whose drizzle weight is 0 is excluded. The spec needs to say that:
   - w in W, W2, the history test's w_new and the LSF_REJ counts is the per-pixel product;
   - COUNT counts samples with d > 0, per channel.
9. **Drizzle weights must be comparable across appends (spec + code).** Siril renormalizes the drizzle weights per channel by `driz->max_weight[c]`, computed per apply run from the reference and the optional flat (applyreg.c `compute_max_drizzle_weights`). By default it also saves them as 8-bit (`drizz_weight_match_bitpix`). That is the §10.4 problem again: per-run normalization makes batches incomparable. LSF must:
   - either use unnormalized weights or freeze the normalization at init;
   - force float32 weight files;
   - record per-frame drizzle facts (pixfrac, kernel, max_weight) in the ledger.
10. **§7 instrument checks under drizzle (spec).** "XPIXSZ, YPIXSZ and FOCALLEN, scaled for drizzle the way Siril does it" is set at init and then checked against incoming frames. Raw CFA frames carry the unscaled pixel size, so the check must apply the same scaling before it compares. Today the comparison would refuse every frame of a 2x file.
11. **'wcs' registration needs a WCS on the Stack, and nothing writes one (spec + code).** §7 says WCS is updated "when the stack is plate-solved again". But the LSF commands are the only writers, and none of them can plate-solve. Plate-solving the primary in the GUI and saving goes through Siril's normal save, which drops the extensions. Today a WCS only exists if the first batch was plate-solved before lsf_init, so that the apply step carried it. M5 needs a way to set the WCS (B5).
12. **REG_XFRM cannot describe a distorted registration (spec, §9).** Siril's astrometric apply can undistort with SIP (`init_disto_data`), and the result is not a 3×3 homography. Proposal:
    - REG_TYPE = 'WCS';
    - REG_XFRM holds the linear part, which is enough for the footprint and expansion shifts;
    - the exact map is defined by the frame's own WCS, already stored verbatim in ORIGHDR, plus the grid WCS.

    State this in §9 and §10.6. An expansion shifts the grid's CRPIX, and the frame's WCS is unaffected.
13. **§3 sizes ignore drizzle (spec).** 2x drizzle is 4× the pixels. A 26 MP OSC sensor at 2x gives 104 MP, about 5.6 GB of history in RAM and about 6.9 GB uncompressed on disk. This is the real driver for band-by-band, more than the 26 MP threshold in the milestone note.
14. **Check 3 and optional HDUs (spec, §8).** Check 3 validates only the required HDUs. Once LSF_REJ is real, it should also check LSF_REJ's dimensions, named planes and lossless compression when the HDU is present.

### A3. Minor, can wait

15. LSF_STARS (§10.5) is not implemented. The writer drops it on rewrite, which the spec allows. Simplify §10.6 to "LSF_STARS is discarded" (drop "shifted the same way").
16. `ledger.orig_baseline` / ORIGHDR diffs (§6, §9) are specified but not implemented: the writer always stores verbatim headers. That is fine because the diffs are conditional, but the text could say MAY.
17. §13 still has the open question on the clip correction for percentile/GESD. It is not touched by M5.
18. §12 test plan: add CFADRIZ parity, band-by-band parity and 'wcs' registration (see B7).

## B. Work items

Order matters: B3 (band-by-band) and B4 (CFADRIZ) both rewrite the core of `lsf_process_batch` and the writer, so they run one after the other, not in parallel.

### B1. Spec update (documentation only, first)

Apply A1.1, A1.4–A1.6 and A2.8–A2.14 to "live stack fits.md", plus the outcome of the C-decisions. Keep LSFVER at 0.4, because nothing here changes the layout of existing files except adding optional ledger columns and an optional HDU.

### B2. LSF_REJ (if C3 = yes)

- `lsf_init -rejmaps`; `lsf_append -rejmaps` enables it later and sets LSFRGEN. The HDU's presence is the fact, and `rejmaps.enabled` mirrors it.
- Keep a per-pixel uint16 low/high counter per channel in `lsf_history` (2·P planes, +4 B/px/chan), fed by both the batch-local and the history test. It saturates and sets LSFRSAT.
- Writer: RICE_1 when compressed, LSFPLn names, WCS copies; read and preserve it on rewrite, and pad it on expansion.
- `lsf_maps -rejlo -rejhi`.
- Fixes the A1.2 drop bug.
- Test: counts match the HISTORY percentages × the number of pixels.

### B3. Band-by-band history I/O

Goal: peak memory = one band of history + the batch's frame buffers, independent of grid size. One code path: "in memory" is just a single band of full height.

- **Read.** Load rows [y0, y1) of each of the 4P (+2P) planes. The compressed tiles are one row each (§5.3), so row-band reads stay cheap. This is a good reason to keep that tiling, and §5.3 should say so.
- **Merge.** Wrap the existing block loop in an outer loop over bands. The frame regions read per block are already row-bounded.
- **Write, pass 1.** Create HDU 0 with its final header, using placeholder LSFHASH/LSFNMIN/LSFNMAX values of fixed width, and leave its data unwritten. Write each band of LSF_ACC, LSF_COUNT (and LSF_REJ) into the temporary file. Accumulate the min/max of MEAN where W > 0 for output normalization.
- **Write, pass 2.** Re-read MEAN/W bands from the temporary file, build HDU 0 bands in file order, hash them incrementally, then fill in the placeholder keys. Without outnorm, HDU 0 can be written during pass 1.
- **Grid expansion** while streaming: output band rows map to input rows offset by pad; pad rows and columns are zero.
- **Spike first:** out-of-order tile writes into a compressed 3-D cube. A band loop writes plane 1 rows a..b, then plane 2 rows a..b, which is not tile order. If CFITSIO cannot do this, write each plane's bands sequentially and re-read the frame regions per plane (slower), or stage the compressed HDUs.
- **Also stream** check 8 (the hash over HDU 0), `-repair`/`-adopt`/`-compress` in lsf_util, and lsf_maps (one output map in RAM at a time).
- **Registration** reads all of HDU 0 as the external reference: 1.25 GB for a 2x OSC grid. Consider loading only the registration layer. This needs a small applyreg change, or a temporary single-layer reference file.
- **Gate:** bit-identical Stack, ACC, COUNT and ledger against the current in-memory build on the existing test sets, with a forced small band height (for example 37 rows), compressed and uncompressed, with and without an expansion. This is the same method as the compression parity test (§12 check 12).

### B4. CFADRIZ

Stage 1 is scale 1 (Bayer drizzle without upsampling), which fits in memory today. Stage 2 is scale > 1, after B3.

- **Init.** `lsf_init -drizzle [-scale=] [-pixfrac=] [-kernel=] [-flat=]` on a calibrated, non-debayered CFA sequence with registration computed. Apply it with `regargs->driz` (Bayer) onto the canvas. Check that the `regargs->canvas` patch scales correctly with `output_scale`. Then LSFDMOD = 'CFADRIZ', `cfa.mode`, `grid.oversample`.
- **Weights.** Read the drizztmp weight files per block, as Siril's stack does (`dstack`). Sample weight = frame weight × d, d = 0 excluded, and per-pixel w_new in the history test (A2.8).
- **Weight scale.** Force float32 weight files and freeze or avoid the `max_weight` renormalization (A2.9). Spike: confirm what `max_weight` depends on (scale, pixfrac, kernel, flat), and pick between freezing it at init and dividing it out.
- **Append.** Register CFA frames against the RGB Stack. If the Stack is at scale s, the star-alignment H already contains s, and the drizzle apply multiplies by s again, so divide before applying. Spike on a 2x file.
- **Ledger.** New optional columns DRZ_PIXF, DRZ_KERN and DRZ_WMAX (per plane). pixfrac and kernel become per-append options with defaults, not settings (A1.5).
- **§7 checks** with drizzle scaling (A2.10).
- **Out of the box:** mono drizzle (cfa.mode NONE, oversample > 1) uses the same path. Include it or refuse it explicitly (C5).
- **Gate:** parity with `seqapplyreg -drizzle` + `stack` on a CFA set (ATR2600C or Nikon data), bit for bit under the §12 check 11 rules, at scale 1 and 2.

### B5. grid.reference = 'wcs' (astrometric registration)

- **Stack WCS:** add `lsf_util -wcs=solved.fits`. It imports the WCS cards from a plate-solved image that must have the Stack's dimensions, and preferably the lsf_export provenance HISTORY with the same LSFUID. It rewrites them into HDU 0 and every extension, adds a HISTORY card and increments LSFGEN. The data is unchanged, so the hash is too. The workflow becomes lsf_export → platesolve → lsf_util -wcs. The alternative, calling the solver from lsf_util, is more work and duplicates platesolve.
- **Append:** `lsf_append -wcsreg` (per append, not stored, A1.5). The frames must be plate-solved (seqplatesolve). Use `compute_Hs_from_astrometry` with the Stack's WCS as the reference. This needs an applyreg patch similar to the canvas one: an external WCS reference instead of a frame of the sequence. The footprint, expansion and canvas logic stay as they are.
- **Distortion:** REG_TYPE = 'WCS' with the linear part in REG_XFRM (A2.12). Decide whether to support undistortion at all in M5, or only linear astrometric registration.
- **Gate:** the same batch registered by stars and by WCS gives transforms within about 0.1 px, and §12 check 10 (expansion keeps WCS and ledger in agreement) passes.

### B6. v0.3 → v0.4 converter (if C2 = keep)

The best place is Python in lsf_poc: lsf_fits.py already reads v0.3, and C does not need the old layout. A C version would be `lsf_util -upgrade`, if Siril users are meant to run it.

- Steps: invert outnorm to get MEAN, pack LSF_ACC, compute LSFHASH, map the settings rows (add norm.* if missing), set LSFVER 0.4, keep the ledger.
- Verify by opening the result with `lsf_util -check`, then appending one batch.

### B7. Tests (§12 additions)

- band-by-band parity (B3 gate);
- CFADRIZ parity at scale 1 and 2;
- a CFADRIZ append with an expansion;
- 'wcs' versus star registration agreement;
- LSF_REJ counts versus the HISTORY percentages;
- the converter round trip;
- a 2x drizzle file in the external-tool checks (already listed in §12).

## C. Decisions needed

1. **In-place append (§10.3):** DECIDED 2026-10-22: dropped. Delete §10.3. Keep LSFBUSY as a reserved keyword (check 2 stays, so check numbering is unchanged).
2. **v0.3 converter:** do you have v0.3 files worth keeping? If not, drop it and delete the §10.11 paragraph (recommended). If yes, write it in Python in lsf_poc.
3. **LSF_REJ:** pull it into M5 (recommended: it is cheap, B3 has to design around it anyway, and it fixes the check 7 bug)?
4. **grid.reference, pixfrac, kernel:** per-append options plus ledger facts rather than stored settings (recommended, per the facts-not-policy rule)?
5. **Mono drizzle:** support it alongside CFADRIZ (nearly free), or refuse it?

Suggested order: B1 → B2 → B3 → B4 stage 1 → B4 stage 2 → B5 → B6. B4 stage 1 could go before B3 if CFADRIZ is wanted sooner, at the cost of touching the merge loop twice.
