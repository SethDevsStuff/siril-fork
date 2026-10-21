# Living Stack Format (LSF) - FITS Profile

Version 0.4 (draft). Serialization profile of the LSF data model for FITS. Companion to the XISF profile (v0.2).

## 1. Scope and Priority

This document defines how an LSF file is stored as a standard FITS file. The LSF data model follows v0.2 and is not repeated here: the weighted accumulator (v0.2 Section 5), canonical grid (Section 7), frame ledger (Section 8), ingestion-time normalization (Section 8.3), operations (Section 9) and sensor handling (Section 10). This profile differs from v0.2 in four ways:

- the container;
- an added accumulator plane, the sum of squared weights (W2);
- a raw-unit copy of the mean (MEAN), stored with the accumulators, which is the authoritative mean;
- the weighting, registration and rejection rules of Section 10.

The v0.2 data model needs the same additions.

**Governing priority: a stock FITS reader must open an LSF file and see the current stack, correctly, without knowing LSF exists.** Wherever LSF convenience conflicts with this, compatibility wins. The words MUST, SHOULD and MAY are used as in RFC 2119.

**Design principle: augment Siril, do not change it.** LSF reuses Siril's existing calibration, registration, star detection, normalization, rejection and weighting code. It adds only the accumulator merge, the history rejection test, the frame-weight scaling of Section 10.4 and the file I/O. LSF operations are exposed only through dedicated LSF commands: lsf\_init, lsf\_append, lsf\_export, lsf\_check and lsf\_maps. Existing commands, stacking paths and save paths are unchanged.

**Current scope.** This is a proof of concept. Only Siril appends to an LSF file. Other tools are expected to read the file, not to keep the full LSF structure when they save it (see Section 15).

## 2. Compatibility Rules

- **C1. HDU 0 is the real Stack.** It is a plain, uncompressed, full-resolution float32 image, with WCS whenever the grid has one. It is never a thumbnail and never an empty placeholder.
- **C2. Standard structures only.** Extensions are standard IMAGE or BINTABLE HDUs. No custom XTENSION types, no non-standard BITPIX, no HIERARCH keywords required for correct operation.
- **C3. Standard keywords stay true.** Every standard keyword in the primary header describes the Stack as it is now, and is updated on every append (Section 7). A keyword whose value can no longer be kept true is removed, not left wrong.
- **C4. The image never depends on an extension.** A reader that discards every extension loses LSF abilities, not the image.
- **C5. Self-identifying and tamper-evident.** LSF keywords use the LSF prefix and stay within 8 characters. Extension names use the LSF\_ prefix. The primary is bound to the extensions by an integrity check (Section 8), so foreign edits are detected rather than silently merged.
- **C6. Locate by name, never by position.** Readers find HDUs by EXTNAME, and planes within a cube by their LSFPLn keywords.
- **C7. Preserve the unknown.** Writers MUST preserve unrecognized HDUs, keywords and table columns when rewriting a file.

## 3. File Layout

| Order | Type | EXTNAME | Contents | v0.2 (XISF) equivalent |
| --- | --- | --- | --- | --- |
| 0 | PRIMARY | none | Stack as displayed: the weighted mean, output-normalized when LSFNORM = T, float32 | Stack |
| 1 | IMAGE | LSF\_ACC | Accumulator cube, float32: MEAN, W, M2 and W2, P planes each | Stack (raw), WeightSum, SumSq, none (W2 is new in v0.3) |
| 2 | IMAGE | LSF\_COUNT | Contribution count, uint16 | ContributionCount |
| optional | IMAGE | LSF\_REJ | Rejection counts, low and high, uint16. Off by default (Section 5.4) | none |
| 3 | BINTABLE | LSF\_SETTINGS | Settings, grid and normalization parameters | Settings, CanonicalGrid, OutputNorm |
| optional | BINTABLE | LSF\_STARS | Cached registration star list (Section 10.5) | none |
| last | BINTABLE | LSF\_LEDGER | One row per contributing frame | FrameLedger |

P is the number of channels: 1 for mono, 3 for OSC (R, G, B).

- **LSF\_ACC** has NAXIS = 3 and NAXIS3 = 4·P. Planes are ordered MEAN, W, M2, W2, each with P channels in R, G, B order. For OSC that is MEAN\_R, MEAN\_G, MEAN\_B, W\_R, ... W2\_B. Each plane is named by an LSFPLn keyword in the LSF\_ACC header, for example LSFPL1 = 'MEAN\_R'; mono files use 'MEAN', 'W', 'M2' and 'W2'. Readers find planes by these names (C6).
  - MEAN is the weighted mean in raw accumulator units, and is the authoritative copy of the mean.
  - W is the sum of weights, M2 the Welford M2 and W2 the sum of squared weights.
  - All four planes are always per channel.
- **Why one cube.** Siril treats a FITS file holding several image HDUs of the same size and BITPIX as a FITS sequence (src/io/fits\_sequence.c, \_find\_hdus). An HDU whose NAXIS3 or BITPIX differs from the first image is skipped. With separate float32 planes, stock Siril opened an LSF file as a 4-frame sequence. NAXIS3 = 4·P (4 or 12) never equals the primary's NAXIS3 (none or 3), and LSF\_COUNT and LSF\_REJ differ in BITPIX. Stock Siril therefore sees a single image, without any change to its reader.
- **LSF\_COUNT** has the same grid as the Stack. Mono files use NAXIS = 2. OSC files use NAXIS = 3 with NAXIS3 = 3. LSF\_COUNT may use NAXIS3 = 1 only when every channel always accepts or rejects the same samples (for example cfa.mode = DEBAYER with rejection disabled). It never may for CFADRIZ, where coverage differs per channel.
- **HDU 0 is derived from MEAN.** It is rebuilt from MEAN on every commit: output-normalized when LSFNORM = T, otherwise clamped to \[0, 1\] (Section 4.1), and 0 wherever W = 0. HDU 0 is load-bearing for compatibility (C1), but it is no longer the only copy of the mean, so a damaged primary can be repaired (Section 8).
- **Ledger placement.** An LSF writer puts LSF\_LEDGER last. Readers MUST NOT require it to be last (C6). When rewriting, a writer places unknown HDUs (C7) after the LSF image and settings HDUs and before the ledger, in their original relative order. Optional LSF\_ HDUs, including the deferred per-frame data of v0.2 Section 14, also go before the ledger.

**Size.** Per pixel per channel: Stack 4 B, LSF\_ACC 16 B, LSF\_COUNT 2 B, total 22 B. With LSF\_REJ, add 4 B, for 26 B. On a 26-megapixel sensor this is about 575 MB mono and 1.72 GB OSC, or 680 MB and 2.03 GB with rejection maps. Lossless compression of the extensions (Section 5.3) typically saves 20 to 35 percent of the extension size.

## 4. Primary HDU

### 4.1 Image

- BITPIX = -32. NAXIS = 2 (mono) or 3 with NAXIS3 = 3 (OSC, planes R, G, B).
- BSCALE and BZERO MUST be absent or 1 and 0 (Siril writes BSCALE = 1 and BZERO = 0).
- The primary MUST carry ROWORDER, as Siril does ('TOP-DOWN' in the sample header). Every extension image MUST use the same row order, because the registration transforms and any WCS depend on it.
- Pixels with no data (weight 0) are written as 0.0, not NaN. NaN breaks statistics and auto-stretch in many viewers, and 0 matches the output-normalization behavior in v0.2 Section 9.
- When output normalization is enabled, the primary is MEAN with the same min-max output normalization Siril applies to stacks, and is in \[0, 1\]. Otherwise it is MEAN clamped to \[0, 1\], as Siril clamps its own 32-bit stacks. With the same frames and settings, the primary then equals the image Siril's stack command produces. MEAN is never clamped, so values above 1 (star cores raised by normalization) stay available to the accumulator and to a later output normalization.
- An OSC primary is already RGB data. BAYERPAT, XBAYROFF and YBAYROFF MUST NOT appear in the primary header, because tools that see BAYERPAT may try to debayer an already-debayered image. The pattern is recorded in LSFBAYR instead. This matches Siril 1.4.4, whose OSC stack header carries none of these keywords.

### 4.2 LSF keywords

| Keyword | Type | Meaning |
| --- | --- | --- |
| LSF | logical | T. Marks the file as LSF |
| LSFVER | string | Profile version, e.g. '0.4' |
| LSFUID | string | Unique id, created at initialization, never changed |
| LSFGEN | integer | Generation counter, incremented on every commit, identical in all LSF HDUs |
| LSFSENS | string | 'MONO' or 'CFA' (permanent, v0.2 Section 10.1) |
| LSFBAYR | string | Bayer pattern (RGGB, GRBG, GBRG, BGGR). CFA only |
| LSFDMOD | string | 'NONE' (mono), 'DEBAYER' or 'CFADRIZ' (permanent, v0.2 Section 10.3) |
| LSFNORM | logical | T if the Stack is output-normalized |
| LSFNMIN, LSFNMAX | float | Min and max of the output-normalization map |
| LSFHASH | string | SHA-256 of the primary data at last commit, lowercase hex (Section 8) |
| LSFBUSY | logical | Present and T only during an in-place update (Section 10.3) |

LSFHASH is computed over the primary data unit as stored: big-endian float32 values in file order, without the trailing padding to a 2880-byte block. Its 64 characters fit in a string card.

LSFNMIN and LSFNMAX MUST be written with enough digits to round-trip exactly: 9 significant digits for float32 values, 17 for float64. Repair (Section 8) uses them to rebuild the primary from MEAN, and adopt (Section 8) uses them to invert the normalization. A writer MUST refuse to adopt a primary flagged LSFNORM = T if either is missing.

## 5. Extension Image HDUs

### 5.1 Headers

Each extension image carries XTENSION = 'IMAGE', BITPIX, NAXIS and axes, PCOUNT = 0, GCOUNT = 1, EXTNAME, EXTVER = 1, LSFUID and LSFGEN. Cubes (LSF\_ACC, and LSF\_REJ or a per-channel LSF\_COUNT) also carry an LSFPLn keyword naming each plane.

Each also carries **its own copy of the WCS keywords** (CTYPEn, CRVALn, CRPIXn, CD matrix, any distortion terms), but only when the grid has a WCS. INHERIT is not used, because readers support it unevenly. For a cube, the WCS keywords describe axes 1 and 2 only.

DATASUM and CHECKSUM are RECOMMENDED on every HDU. Siril's LSF commands always write them.

### 5.2 Count plane

LSF\_COUNT uses BITPIX = 16 with BSCALE = 1 and BZERO = 32768, the standard unsigned 16-bit convention. The count **saturates at 65535 and never wraps**. A writer that saturates any pixel sets LSFCSAT = T in the LSF\_COUNT header.

Under reliability weights, saturation does not affect n\_eff, which is computed from W2 (Section 10.8). Under noise weighting the count enters the variance estimate as N. There, a saturated pixel's estimate is approximate, with a relative error below 1/65535. Otherwise a saturated count only makes the coverage information at that pixel approximate.

### 5.3 Compression

- **Default.** Compression is optional and off by default. It is enabled by an option on the LSF commands, as with the XISF reference writer.
- **What may be compressed.** Extension images MAY be tile-compressed (ZIMAGE convention). Siril's FITS compression support is the reference for which methods and parameters are available, and the lossless requirement below still applies.
- **Lossless only.**
  - LSF\_ACC uses GZIP\_2 (byte-shuffled) with float quantization disabled. RICE, HCOMPRESS and any quantized float mode are prohibited for it.
  - LSF\_COUNT and LSF\_REJ may use RICE\_1 or GZIP\_2.
- **Primary.** The primary HDU MUST NOT be tile-compressed, and MUST NOT be replaced by an empty primary followed by a compressed image extension. That pattern breaks C1. Siril's existing compressed save writes exactly this pattern, so the LSF writer cannot use it for HDU 0.
- **Whole-file gzip** (.fits.gz) is allowed for transport only. A file must be decompressed before appending.
- **Updates.** Compressed HDUs cannot be updated band by band in place. A writer using compression rewrites them on append (Section 10.2). Expect roughly the 20 to 35 percent saving described in v0.2 Section 6.1.
- **XISF modes.** The XISF-specific modes (zlib+sh, zlib-mt, subblocks) do not carry over. CFITSIO tile compression replaces them.

### 5.4 Rejection planes (optional)

LSF\_REJ is off by default. It is created only when the user asks for it, by an option on lsf\_init or lsf\_append. It records, per pixel and channel, how many samples were rejected low and how many high.

- It uses BITPIX = 16 with BSCALE = 1 and BZERO = 32768.
- NAXIS3 = 2·P, ordered REJLO (R, G, B) then REJHI (R, G, B), with planes named by LSFPLn ('REJLO\_R' and so on; 'REJLO' and 'REJHI' for mono).
- Low and high are always stored separately, never merged. The two directions mean different things: high rejections come from satellites, aircraft, hot pixels and cosmic rays, low rejections from cloud, cold pixels and changing dust shadows.
- It counts rejections from both the batch-local test and the history test (Section 10.7).
- Counts saturate at 65535 and never wrap. A writer that saturates any pixel sets LSFRSAT = T in the LSF\_REJ header.
- When rejection maps are enabled on a later append rather than at initialization, they count only from that commit on. LSFRGEN in the LSF\_REJ header records the first generation counted.
- The setting rejmaps.enabled mirrors whether the HDU exists.

## 6. LSF\_SETTINGS

A BINTABLE with one row per setting: KEY (32A), TYPE (8A: STR, INT, FLT, BOOL or JSON) and VALUE (1PA, variable-length string). A table is used instead of header keywords because settings include nested structures (rejection parameters, alignment options), and new settings should not need a schema change.

Required keys:

- sensor.class, sensor.bayer, sensor.binning, cfa.mode
- grid.oversample, grid.pixfrac, grid.kernel
  - The grid size and WCS are the primary's NAXIS and WCS keywords, and are not duplicated.
- grid.reference ('stack' or 'wcs', Section 10.5), grid.origin\_x, grid.origin\_y (total shift from grid expansions, Section 10.6), grid.expand, grid.expand\_step
- align.transform, align.interpolation, align.min\_pairs, align.findstar (JSON: the star-detection parameters used for registration)
- reject.method, reject.params, reject.hist\_low, reject.hist\_high, reject.hist\_margin (Section 10.7)
- weight.scheme, plus the frozen weight constants for that scheme (Section 10.4):
  - noise: weight.noise\_norm
  - wfwhm: weight.fwhm\_min, weight.fwhm\_max, weight.fwhm\_norm
  - nbstars: weight.stars\_min, weight.stars\_max, weight.stars\_norm
- norm.method
- outnorm.enabled, outnorm.min, outnorm.max
- rejmaps.enabled
- filter (mono and OSC files), created\_utc, writer (name and version)

Conditional key: ledger.orig\_baseline, the baseline header that ORIGHDR diffs are taken against (Section 9). Required when any ledger row stores a diff.

Settings are recalled on append (v0.2 Section 6). Some values are mirrored in the primary header: sensor.class and LSFSENS, sensor.bayer and LSFBAYR, cfa.mode and LSFDMOD, outnorm.\* and LSFNORM, LSFNMIN, LSFNMAX, and rejmaps.enabled and the presence of LSF\_REJ. These MUST agree. A disagreement fails validation (Section 8).

## 7. Standard Keywords

The primary header follows what Siril's own stacks write (STACKCNT, LIVETIME, EXPSTART and EXPEND in Siril's stacking code). Under C3, each keyword below falls into one of three groups.

**Updated on every append:**

| Keyword | Value |
| --- | --- |
| STACKCNT | Number of frames in the ledger |
| LIVETIME | Sum of EXPTIME over the frames in the ledger, seconds |
| EXPSTART | Julian date of the earliest frame start, as Siril computes it for stacks |
| EXPEND | Julian date of the latest frame end (start plus exposure), as Siril computes it for stacks |
| DATE-OBS | Earliest frame start, UTC |
| DATE | Time the file was written, UTC |
| PROGRAM | Writer name and version |
| HISTORY | One card per append, giving the generation, batch number and frame count, followed by one card per channel with the percentages of samples rejected low and high by the batch-local and the history tests. Repair and Adopt (Section 8) add one card each |
| WCS keywords | Updated on grid expansion (CRPIX, Section 10.6) or when the stack is plate-solved again |
| MIPS-HI, MIPS-LO | Removed on append, because display cutoffs saved for the previous stack no longer fit |

**Set at initialization, then checked on append.** A mismatch fails the compatibility check and the frame is not merged:

- XBINNING, YBINNING
- XPIXSZ, YPIXSZ and FOCALLEN, scaled for drizzle the way Siril does it
- INSTRUME, FILTER

**Copied from the first frame.** If a later frame differs, the keyword is removed from the primary, and the ledger keeps the per-frame values:

- EXPTIME (the sub-exposure length). LIVETIME stays exact when it is removed.
- GAIN, OFFSET, CCD-TEMP, EGAIN
- TELESCOP, OBJECT

**Never in the primary:** BAYERPAT, XBAYROFF, YBAYROFF (Section 4.1).

## 8. Integrity Binding

The primary is a view of MEAN. A foreign tool that edits it while preserving the extensions would produce a file whose image no longer matches its accumulators. FITS has no way to prevent this. The profile makes it **detectable** and **repairable**.

**Binding:** LSFUID and LSFGEN in every LSF HDU, LSFHASH in the primary, and DATASUM and CHECKSUM on every HDU.

LSFHASH replaces the LSFDSUM of v0.3. DATASUM is a ones-complement sum of 32-bit words, and its result does not depend on word order. A tool that flips, mirrors or rotates the primary by 180 degrees leaves DATASUM unchanged, and a DATASUM comparison would accept the file. SHA-256 depends on position. DATASUM and CHECKSUM stay as the standard transport checks.

**What the hash does.** It detects; it cannot prevent. It cannot stop another program overwriting the file. It makes Siril's LSF commands refuse to merge into a primary that no longer matches MEAN. Within Siril, the LSF commands are the only path that writes an LSF file, and every write goes through the commit protocol (Section 10.2).

**Pre-merge validation.** A writer MUST check, in order:

1. LSF = T and LSFVER is supported (Section 11).
2. LSFBUSY is absent or F.
3. Every required HDU is present (PRIMARY, LSF\_ACC, LSF\_COUNT, LSF\_SETTINGS, LSF\_LEDGER, found by EXTNAME) with consistent dimensions, and LSF\_ACC holds the expected named planes.
4. LSFUID and LSFGEN are identical across all LSF HDUs.
5. DATASUM and CHECKSUM match on every extension where they are present.
6. Every ledger GEN is at most LSFGEN. A generation with no ledger rows is one that merged no frame (Repair or Adopt), and is recorded in a HISTORY card.
7. Mirrored settings agree.
8. The SHA-256 recomputed from the primary data equals LSFHASH.

If any of checks 1 to 7 fails, the writer refuses to merge, reports which check failed, and leaves the file untouched. The file stays readable as an image (C1). There is no silent repair. If LSF\_ACC itself is damaged (check 3 or 5), the accumulator cannot be trusted, and the user restores a backup copy.

If only check 8 fails, the accumulators are intact but the primary was changed outside Siril. The writer stops, reports that the Stack was modified outside Siril, and offers two explicit actions. In Siril they are options of lsf\_check, which otherwise only validates the file and summarizes it:

- **Repair** (the default, lsf\_check -repair). Rebuild HDU 0 from MEAN (Section 3). This undoes the outside edit: the rebuilt primary is bit-identical to the one written at the last commit.
- **Adopt** (lsf\_check -adopt). Accept the edited primary as the new mean. A pixel whose primary value is still the one its MEAN gives was not edited and keeps MEAN exactly, which also keeps values the primary clamps (Section 4.1). An edited pixel gets MEAN from HDU 0, inverting output normalization with LSFNMIN and LSFNMAX when LSFNORM = T. Pixels with W = 0 stay 0. M2 is not adjusted, so adopt is meant for small cosmetic fixes the user wants to keep. Adopt cannot work if the dimensions have changed. The primary is then rebuilt from the new MEAN, so with output normalization an edit beyond the old range rescales the whole Stack.

Either action records a HISTORY card, increments LSFGEN, and recomputes LSFHASH and the checksums.

**Expected behavior of foreign tools:**

| What the other tool does | Result |
| --- | --- |
| Reads only HDU 0 | Correct current stack |
| Lists all image HDUs | Sees LSF\_ACC as a 4·P-plane cube and LSF\_COUNT. The LSF\_ names show they are not frames. Stock Siril skips them (Section 3) |
| Saves only the primary | A plain FITS file. Extensions are gone, so check 3 fails. The image is intact |
| Strips unknown keywords | The file stops identifying as LSF and becomes a plain FITS |
| Re-saves with extensions but a changed primary | Check 8 fails. Repair restores the primary from MEAN, or the user adopts the edit |
| Flips, mirrors or rotates the primary | Check 8 fails. A DATASUM comparison would not have detected it |
| Reorders HDUs | No effect (C6) |

These are design expectations, not observed behavior, except for the Siril behavior in the second row, which was confirmed against the Siril source. See Section 12.

## 9. LSF\_LEDGER

A BINTABLE, written as the last HDU (Section 3), with one row per contributing frame, carrying the fields of v0.2 Section 8. Per-plane columns have repeat count P. Their TFORM is written below as nD or nE, where n = P: 1D for mono, 3D for OSC.

| Column | Format | Content |
| --- | --- | --- |
| FRAME\_ID | 32A | Unique id |
| BATCH | J | Batch number |
| GEN | J | LSFGEN of the commit that added the frame |
| ADDED\_UTC | 24A | Ingestion time |
| FILENAME | 1PA | Original file name |
| DATE\_OBS | 24A | Exposure start, UTC, ISO 8601 |
| EXPTIME | D | Exposure, seconds |
| GAIN, OFFSET | J | Sensor settings (TNULL when unknown) |
| XBINNING, YBINNING | J | Binning |
| CCD\_TEMP | E | Sensor temperature, C |
| FILTER | 16A | Filter |
| ORIGHDR | 1PA | Verbatim source FITS header, or a diff against ledger.orig\_baseline (Section 6) |
| REG\_TYPE | 8A | Transform type |
| REG\_XFRM | 9D | 3x3 transform onto the current grid, row-major, TDIM = (3,3) |
| ALIGN | 1PA | Alignment parameters used (JSON) |
| NORM\_SC, NORM\_OFF | nD | Normalization scale and offset per plane |
| NORM\_MTH | (12·P)A, TDIM = (12,P) | linear, offset\_only, failed or skipped, per plane (12A mono, 36A OSC) |
| WEIGHT | nD | Weight applied per plane, after the scaling of Section 10.4 |
| FWHM, WFWHM, ROUND, REGQ, REGBG | E | Quality metrics |
| NSTARS | J | Star count |
| BKG\_MED, NOISE, NOISE\_N | nE | Per-plane background and noise |
| STAT\_MEAN, STAT\_SIG, STAT\_MAD, STAT\_SQBW, STAT\_LOC, STAT\_SCL | nE | Per-plane Siril image statistics |

Rules:

- Readers locate columns by TTYPE. Every quality column is optional. Unknown values are NaN (floats), TNULL (integers) or empty (strings). Files written by earlier versions simply lack columns.
- Unknown columns MUST be preserved on rewrite (C7).
- String columns hold printable ASCII only (32 to 126), as FITS requires. ORIGHDR is stored the way FITS stores a header: 80-character cards one after another, padded with spaces, with no line breaks. A reader splits it every 80 characters. Siril keeps headers as text with one card per line, so the writer removes the line breaks.
- The ledger is append-only in normal operation. Because the writer puts it last, the writer can rewrite it in full on every append without moving any other HDU. It holds scalars and short text, so this stays cheap.
- REG\_XFRM always maps onto the grid as it is now. A grid expansion rewrites the existing rows (Section 10.6).

## 10. Operations

### 10.1 Initialize

lsf\_init does the following:

1. Builds the canonical grid from the first batch's reference frame, using Siril's existing registration.
2. Merges the first batch (v0.2 Section 9).
3. Freezes the weight constants (Section 10.4).
4. Writes every HDU. LSF\_REJ is written only if the user asked for rejection maps.
5. Sets a new LSFUID and LSFGEN = 1, and writes LSFHASH and the checksums.

If the first batch has fewer than 25 frames, the writer SHOULD warn the user before creating the file. Per-pixel scatter estimated from so few frames is unreliable, and any rejection decision made on it stays in the stack. The warning does not block creation. The profile defines no minimum stack size for rejection, so rejection behaves the same at every count. The threshold of 25 is a default.

### 10.2 Append (preferred commit protocol)

1. Run the Section 8 validation.
2. Process the batch with Siril's existing code:
   - the v0.2 compatibility check, which includes the Section 7 keywords that are set at initialization;
   - calibration;
   - registration against the grid reference (Section 10.5);
   - normalization (v0.2 Section 8.3).
3. Grow the grid first if needed (Section 10.6).
4. Compute the frame weights (Section 10.4) and run rejection (Section 10.7).
5. Update the accumulator, in memory or band by band.
   - The batch is merged into history with the weighted pairwise update. With δ = MEAN\_b − MEAN\_a and W = W\_a + W\_b:
     - MEAN = MEAN\_a + δ·W\_b / W
     - M2 = M2\_a + M2\_b + δ²·W\_a·W\_b / W
     - W2 = W2\_a + W2\_b
     - COUNT = COUNT\_a + COUNT\_b, saturating
   - LSF\_REJ is updated if present.
6. Rebuild HDU 0 from MEAN, update the Section 7 keywords and append the ledger rows. Increment LSFGEN in every LSF HDU, and recompute LSFHASH and the checksums.
7. Write the complete new file to a temporary name in the same directory.
8. Flush to disk, then atomically replace the original (Section 10.11). A writer SHOULD offer to keep the previous file as a backup.

### 10.3 Append (in-place alternative, uncompressed files only)

1. Set LSFBUSY = T and flush the header.
2. Update the planes, the ledger and the keywords.
3. Recompute LSFHASH, and increment LSFGEN in every LSF HDU.
4. Clear LSFBUSY.

A reader that finds LSFBUSY = T knows an update was interrupted. The Stack is still viewable, but merging is refused until the user restores a backup. A writer using in-place updates SHOULD keep a backup.

### 10.4 Frame weights

Siril's existing weighting functions normalize the weights of one stacking run so that their mean is 1. The wFWHM and star-count weights are also scaled between that run's minimum and maximum, so the worst frame of a run gets weight 0. If these were applied to each batch independently, every batch would carry the same weight per frame whatever its quality, and weights would not be comparable across appends.

LSF therefore reuses Siril's per-frame inputs and Siril's formulas, but with constants frozen from the first batch:

- **noise:** the raw weight is 1 / (pscale² · bgnoise²), using Siril's normalization scale and background noise for the channel. weight.noise\_norm is the mean raw weight of the first batch, and every later weight is the raw weight divided by it. No clamping is needed.
- **wfwhm:** Siril's formula, with weight.fwhm\_min and weight.fwhm\_max taken from the first batch, and weight.fwhm\_norm the first batch's mean weight. A later frame's wFWHM is clamped to \[fwhm\_min, fwhm\_max\] before the formula is applied.
- **nbstars:** Siril's formula, with weight.stars\_min, weight.stars\_max and weight.stars\_norm from the first batch. Star counts are clamped the same way.
- **none:** every weight is 1.

A frame whose weight is 0 after clamping contributes nothing. It is not merged, and it is reported to the user. These rules live in an LSF helper alongside Siril's code; Siril's own weighting functions are not changed. Siril's nbstack weighting is not supported in this version.

### 10.5 Registration and the grid reference

The canonical grid is the Stack's own pixel grid. At initialization it is defined by the first batch's reference frame, scaled by grid.oversample. That frame is not stored, so later batches are registered against **HDU 0 itself**. They use Siril's existing star detection, with the parameters stored in align.findstar, and Siril's existing global registration.

The Stack is the highest-SNR reference available and already defines the grid. The writer records grid.reference = 'stack'. When the Stack has a WCS, grid.reference MAY be 'wcs', and new frames are registered with Siril's astrometric registration instead.

The optional LSF\_STARS BINTABLE caches the star list detected on HDU 0, typically under 1 MB, to make appends faster and repeatable. Its columns are X, Y, FLUX and FWHM, in grid pixels. Its header carries LSFGEN. A reader uses the cache only when that LSFGEN equals the file's; otherwise it detects stars again.

### 10.6 Grid expansion

Framing drifts between sessions, so the grid must be able to grow. When grid.expand = T (the default) and accepted frames would cover pixels outside the current grid, the writer grows the grid before merging:

- Growth is in multiples of grid.expand\_step pixels (default 64), so small drifts do not rewrite the file on every append.
- Every image plane is padded with zeros. W = 0 means no data, so padding changes no existing pixel.
- Rows or columns inserted before the existing data, in storage order, are added to CRPIX1 and CRPIX2 in every HDU that carries a WCS, and to grid.origin\_x and grid.origin\_y.
- REG\_XFRM in every existing ledger row is rewritten to include the shift, so every transform keeps mapping onto the current grid.
- LSF\_STARS is shifted the same way, or discarded.

Expansion never shrinks the grid. Cropping is a separate operation, out of scope here.

### 10.7 Rejection

Rejection runs in two stages, in this order:

1. **Batch-local rejection.** Siril's existing rejection is run on the batch, unchanged, with reject.method and reject.params. Siril stops rejecting once 4 or fewer samples remain at a pixel, so a batch of 4 frames or fewer gets no batch-local rejection, and only the history test applies.
2. **History test.** This is the only rejection step that LSF adds. Each surviving sample x of a new frame is compared with the accumulated history. Let σ\_pred be the predicted spread of a new sample about MEAN (Section 10.8). The sample is tested with

   z = (x − MEAN) / σ\_hist, where σ\_hist = reject.hist\_margin · max(c · σ\_pred, σ\_frame).

   It is rejected if z < −reject.hist\_low or z > reject.hist\_high.

A sample is rejected if it fails either stage. The history thresholds are separate settings from the batch-local ones. By default they are equal to the batch sigma thresholds when the batch method has them, and 3.0 otherwise.

**Guarding against a shrinking spread.** M2 only accumulates accepted samples, so the measured spread is biased low. Uncorrected, the history test would reject more on every append, which biases the spread lower still. σ\_hist uses three guards:

- **c**, the clipping correction, is computed from the thresholds alone. For a standard normal clipped to \[−kl, kh\], the remaining variance is

  v = 1 − (kl·φ(kl) + kh·φ(kh)) / Z − ((φ(kl) − φ(kh)) / Z)², where Z = Φ(kh) − Φ(−kl),

  and c = 1 / √v. For each side, k is the tighter of the history threshold and the batch sigma threshold (or the history threshold alone when the batch method has no sigma threshold). For example, c ≈ 1.014 at ±3σ and c ≈ 1.048 at ±2.5σ.
- **σ\_frame** is the incoming frame's noise: Siril's background noise for the channel, multiplied by the normalization scale. It is measured before any rejection, so clipping cannot bias it low. It is a floor, not an estimate.
- **reject.hist\_margin** (default 1.0, at least 1) is a user-tunable widening of the test.

All three are applied when the test runs. Stored M2 stays the exact sum, so the policy can change later without rewriting any file. They need no stored plane: σ\_frame and the frame's weight are already in the ledger (NOISE, WEIGHT), and c comes from the settings.

Where σ\_pred is undefined (fewer than two effective samples), σ\_hist uses σ\_frame alone.

### 10.8 Derived statistics

Siril's existing stacking never computes a weighted per-pixel variance. Its rejection uses unweighted statistics of the pixel stack, and the weights enter only the final mean. The statistics below serve only the history test (Section 10.7) and the exported noise and confidence maps.

With W, M2 and W2 the LSF\_ACC planes and N the LSF\_COUNT value at a pixel and channel, the estimator depends on weight.scheme.

**Reliability weights (none, wfwhm, nbstars).** Here the weights express preference, not noise, and the samples are assumed to share a common variance:

- Effective sample size: n\_eff = W^2 / W2. This replaces the v0.2 noise-map approximation.
- Unbiased variance of one sample: s^2 = M2 / (W − W2 / W). Undefined with fewer than two contributing samples.
- Variance of the mean: s^2 · W2 / W^2.
- Predicted variance of a new sample about the mean: σ\_pred² = s^2 · (1 + W2 / W^2).

**Precision weights (noise).** Here each weight is proportional to 1 / σ² of its frame, and the weights carry the noise model:

- Unit-weight variance: σ0^2 = M2 / (N − 1). Undefined when N < 2.
- Variance of the mean: σ0^2 / W.
- Predicted variance of a new sample with weight w\_new: σ\_pred² = σ0^2 · (1 / w\_new + 1 / W).
- n\_eff = W^2 / W2 is still reported as an effective-integration map.

The two sets agree when all weights are equal. The weight constants are frozen (Section 10.4), so the weights keep the same scale across appends and these estimators stay valid as the stack grows.

### 10.9 Exports

- **Flattened master** (lsf\_export). Copy HDU 0 to a new file and remove the LSF keywords. The result is a plain FITS file. If output normalization is enabled, the primary is already the normalized image. The export keeps a HISTORY card naming the source LSFUID and LSFGEN as provenance.
- **Noise and confidence maps** (lsf\_maps). Derived from the accumulator with the statistics of Section 10.8, never stored, and saved as ordinary FITS images with the primary's header (LSF keywords removed) and a HISTORY card naming the source LSFUID and LSFGEN:
  - **neff**: the effective number of frames, W^2 / W2;
  - **noise**: the standard deviation of the mean;
  - **sigma**: the standard deviation of one frame (of a unit-weight frame under noise weighting);
  - **snr**: MEAN divided by the standard deviation of the mean.

  noise and sigma are in the units of the Stack: divided by LSFNMAX - LSFNMIN when LSFNORM = T. A pixel where a value is undefined (no data, or too few samples) is 0. If LSF\_REJ is present, the low and high rejection maps may be exported the same way.

### 10.10 Undo

The format keeps no per-batch history planes, because they would grow the file with every append. To be able to undo an append, copy the file before appending. To drop a bad night, keep the copy made before it and continue appending to that copy.

### 10.11 Implementation notes

- **Atomic replace on Windows.** A plain rename fails when the target exists, so it needs MoveFileEx with MOVEFILE\_REPLACE\_EXISTING. It still fails while another program, such as a viewer, holds the file open. A file that was just written can also be held for a moment by a virus scanner, so a sharing or access error is retried for a short time first. If it persists, the writer must report this and keep the temporary file, not lose the commit.
- **Long ledger strings and CFITSIO.** CFITSIO 4.7.0 truncates a variable-length string (1PA) that crosses two 2880-byte block boundaries when it is read as text (fits\_read\_col with TSTRING). It returns only the part before the second boundary. ORIGHDR is several kilobytes per row, so some rows are affected in every file, including files written by other tools. A reader that needs the full header reads the column as bytes (TBYTE), using the length from the descriptor (fits\_read\_descript). fitsverify reads these strings the same way and reports "contains non-ASCII text" for the affected rows, because it scans stale bytes in its buffer. That error is a false positive when the stored bytes are printable ASCII. lsf\_compare.py checks the stored bytes directly before it accepts this error.
- **Converting v0.3 files.** A v0.3 file (separate LSF\_WEIGHT, LSF\_M2 and LSF\_W2 HDUs, LSFDSUM) can be converted to v0.4 without loss. MEAN is obtained by inverting the output normalization of the primary, the planes are packed into LSF\_ACC, and LSFHASH is computed.

## 11. Conformance

- **Any FITS reader.** No requirements. It benefits from C1 to C4.
- **LSF reader.** Locates HDUs by name and cube planes by LSFPLn, validates per Section 8, reads the ledger, supports GZIP\_1, GZIP\_2 and RICE\_1 tile compression, and tolerates missing optional keywords, columns and HDUs (LSF\_REJ, LSF\_STARS).
- **LSF writer.** Everything an LSF reader does, plus the pre-merge validation, the commit protocol, preservation of unknown content, and honest maintenance of the primary's standard keywords (Section 7). In Siril, writing is done only by the LSF commands (Section 1).
- **Versioning.** While the major version is 0, a minor version change may be incompatible. A reader accepts only the minor versions it knows, unless a later version of this profile declares itself compatible with an earlier one. From version 1.0 on, a reader accepts a higher minor version, and a new major version may need new parsers. C1 to C4 MUST hold in every version.
- **Files without W2.** A file from before W2 existed can be read and exported, but appending to it is refused, because W2 cannot be rebuilt once samples have been rejected.

## 12. Compatibility Test Plan

Nothing in this profile has been tested against real software yet, except the Siril sequence behavior noted in Section 3, which was confirmed from the Siril source. These are the claims to verify before the profile is frozen.

**Tools:** Siril, PixInsight, DS9, AstroImageJ, ASTAP, astropy or fitsio, fitsverify.

**Test files:**

- mono and OSC;
- uncompressed and compressed extensions;
- with and without LSF\_REJ;
- a 2x drizzle file, to check XPIXSZ against the WCS;
- a file after a grid expansion.

**Checks per tool:**

1. Opens the primary with correct dimensions, values and WCS where one is present.
2. Does not load extension planes as additional frames or channels. For Siril: the file opens as a single image, not a FITS sequence, through both the GUI open and the load command.
3. Plate solving or scale detection agrees with the WCS where one is present, including after a grid expansion.
4. Behavior on Save or Save As: what survives, and whether the file crashes or corrupts.
5. Tolerates compressed extensions after the primary.
6. Tolerates the ledger table's variable-length columns (the heap) without failing or warning. This is a possible compatibility problem with the single-table ledger design.

**LSF behavior checks (Siril):**

7. **Flip detection.** Flip the primary with an outside tool while keeping its old DATASUM. Check 8 must fail, and repair must restore the original.
8. **Shrinking spread.** Run 50 simulated appends from a fixed noise model, with and without the Section 10.7 guards. The rejected fraction per append must stay flat with the guards.
9. **Weight consistency.** Append one good and one poor batch. The total weight each batch adds must reflect its quality, rather than being equal per frame.
10. **Grid expansion.** Append a batch shifted beyond the grid edge. Existing pixels must be unchanged, and the WCS and the ledger transforms must still agree.
11. **Parity with Siril's stack.** lsf\_init over a set of frames must give the same image as Siril's stack command with the same frames, reference and settings, bit for bit, when both normalize with the same frame statistics. Siril caches those in the .seq file with 6 significant digits: the first command on a fresh sequence uses full-precision values and later ones the cached values, so the cache is filled first (seqstat full). Known exceptions: wfwhm and nbstars weighting (Siril keeps a weight-0 frame in the rejection set, LSF drops it), and sigmedian rejection without weights (Siril averages the medians that replaced outliers, LSF counts them as rejected).

## 13. Open Questions

1. **Extra image HDUs (watch item).** Stock Siril was confirmed to load same-sized image HDUs as a FITS sequence. The v0.4 layout avoids this without changing Siril (Section 3). Keep watching for any other target tool that treats the extension images as extra frames, channels or layers, and cover it in the Section 12 tests.
2. **Rejection model.** The v0.2 review asked three questions. This profile answers them in Section 10.7:
   - Order: batch-local rejection runs first, and the history test runs on the samples that survive it. Failing either test rejects the sample.
   - Thresholds: the history test has its own settings.
   - Weighting: the estimator depends on the weighting scheme (Section 10.8).

   Still open:
   - the clipping correction c for batch methods without a sigma threshold (percentile, GESD), which currently fall back to the history thresholds alone;
   - whether grid expansion needs a size cap.

   The shared v0.2 data model needs the same answers.

## 14. Future Improvements

- **Linking multi-filter files.** Filters stay one file per filter. A linking keyword could let related files find each other.
- **Duplicate-frame detection.** Define FRAME\_ID as a content hash of the source data and DATE-OBS, so that appending the same frame twice is refused. Today's stacking has no such guard either.
- **Exposure map and mosaics.** An optional float32 plane holding the per-pixel sum of exposure time over accepted samples (+4 B per pixel per channel). It gives the true integration time when sub lengths differ or for mosaics.
- **Export to FLIS.** FLIS is still in development. The intent is that the Stack becomes the base layer of a FLIS file, with LSFUID recorded in the layer METADATA as provenance, and that the LSF HDUs are not embedded in FLIS.
- **float64 accumulators.** LSF\_ACC with BITPIX = −64, for very deep stacks (+16 B per pixel per channel).
- **nbstack weighting**, for appending stacks rather than single frames.

## 15. Known Limitations

- **Append is Siril-only.** Other tools can open the file and use the primary image, but they are not expected to preserve the extension HDUs, the LSF keywords or the ledger when they save.
- **A foreign save can end the stack's life.** A tool that re-saves the file may drop the extensions or strip the LSF keywords. The result is still a valid image, but the accumulator can no longer be merged into, and continuing the stack means starting again from the original frames. A tool that keeps the extensions but changes the primary is detected and repairable (Section 8).
- **Process exports, not the LSF file.** Edits made to the primary are either undone by repair, or adopted without adjusting M2.
- **Fixed sub-exposure length assumed.** Per-pixel integration time is derived from the count, which is only correct when every frame has the same exposure. Mixed exposures are allowed and LIVETIME stays exact, but per-pixel integration time needs the exposure map (Section 14).
- **No built-in undo.** Undo is done by keeping copies (Section 10.10).
- **Compatibility is unverified.** The behavior of other tools described in this draft is a design expectation until the Section 12 tests are run.
