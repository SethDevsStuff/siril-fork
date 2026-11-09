# Living Stack Format (LSF) - FITS Profile

Version 0.4 (draft). Serialization profile of the LSF data model for FITS. 

## 1. Scope and Priority

This document defines the Living Stack Format (LSF) and how an LSF file is stored as a standard FITS file. An LSF file is a stack that can keep growing: new batches of frames are merged into it without re-reading the frames already in it, and its primary image is always the current stack.

### 1.1 Data model

- **Accumulator.** Per pixel and channel the file keeps the weighted mean (MEAN), the sum of weights (W), the sum of squared deviations from the mean (M2, Welford's incremental algorithm, which avoids the cancellation of raw sums of squares on a bright sky), the sum of squared weights (W2) and the number of accepted samples (COUNT). Once the frames' exposures differ, it also keeps the sum of their exposures (EXP, the exposure map). A batch is merged by updating these planes (Section 10.2); the frames themselves are not kept. Its size depends on the grid, not on the number of frames (Section 3).
- **What it cannot do.** Rejection against the history is a single-pass test of each new sample against the accumulated mean and scatter (Section 10.6). Iterative clipping over the whole history would need the individual frame values, which are not kept. A bad frame that passes per-frame selection and rejection stays in the stack for good, so cosmetic correction and the rejection of clouded, trailed or badly guided frames belong before the append.
- **Canonical grid.** The accumulator lives on a fixed pixel grid with its own WCS, defined by the first batch (Section 10.1). Every later frame is registered and resampled onto it (Section 10.4). The pixel scale and orientation are permanent; the grid can only grow (Section 10.5). Changing the scale or drizzle factor later needs the original frames.
- **Frame ledger.** Every frame that ever contributed has one row: identity, acquisition data, its transform onto the grid, its normalization coefficients, its weight and its quality measures (Section 9). New frames can be graded against this history, and the weights it would have had under another scheme can be recomputed.
- **Normalization.** Each frame is normalized with Siril's own formulas against a baseline measured on the reference frame of the first batch and frozen (norm.baseline, Section 6), so that frames from different sessions stay comparable. Normalization is relative: it keeps the file internally consistent, it is not absolute photometry, and it corrects only what one scale and offset per channel can express. Gradients are removed per frame before the append.
- **Sensors.** A file holds one sensor class, fixed at creation: mono, or CFA with its Bayer pattern. A mono file holds one filter; a multi-filter target is one file per filter. An OSC file keeps every accumulator plane per channel (R, G, B), because the channels differ in noise and background, and is normalized per channel. All channels of a frame share one geometric transform, so atmospheric dispersion cannot be corrected per channel. The demosaic method, debayer before registration or CFA-aware drizzle (cfa.mode), is also fixed at creation; the two are never mixed in one file.

**Governing priority: a stock FITS reader must open an LSF file and see the current stack, correctly, without knowing LSF exists.** Wherever LSF convenience conflicts with this, compatibility wins. The words MUST, SHOULD and MAY are used as in RFC 2119.

**Design principle: augment Siril, do not change it.** LSF reuses Siril's existing calibration, registration, star detection, normalization, rejection and weighting code. It adds only the accumulator merge, the history rejection test, the frame-weight scaling of Section 10.3 and the file I/O. LSF operations are exposed only through dedicated LSF commands: lsf\_init, lsf\_append, lsf\_export, lsf\_util and lsf\_maps. Existing commands, stacking paths and save paths are unchanged.

**Current scope.** This is a proof of concept. Only Siril appends to an LSF file. Other tools are expected to read the file, not to keep the full LSF structure when they save it (see Section 13).

## 2. Compatibility Rules

- **C1. HDU 0 is the real Stack.** It is a plain, uncompressed, full-resolution float32 image, with WCS whenever the grid has one. It is never a thumbnail and never an empty placeholder.
- **C2. Standard structures only.** Extensions are standard IMAGE or BINTABLE HDUs. No custom XTENSION types, no non-standard BITPIX, no HIERARCH keywords required for correct operation.
- **C3. Standard keywords stay true.** Every standard keyword in the primary header describes the Stack as it is now, and is updated on every append (Section 7). A keyword whose value can no longer be kept true is removed, not left wrong.
- **C4. The image never depends on an extension.** A reader that discards every extension loses LSF abilities, not the image.
- **C5. Self-identifying and tamper-evident.** LSF keywords use the LSF prefix and stay within 8 characters. Extension names use the LSF\_ prefix. The primary is bound to the extensions by an integrity check (Section 8), so foreign edits are detected rather than silently merged.
- **C6. Locate by name, never by position.** Readers find HDUs by EXTNAME, and planes within a cube by their LSFPLn keywords.
- **C7. Preserve the unknown.** Writers MUST preserve unrecognized HDUs, keywords and table columns when rewriting a file.

## 3. File Layout

| Order | Type | EXTNAME | Contents |
| --- | --- | --- | --- |
| 0 | PRIMARY | none | Stack as displayed: the weighted mean, output-normalized when LSFNORM = T, float32 |
| 1 | IMAGE | LSF\_ACC | Accumulator cube, float32: MEAN, W, M2 and W2, and the optional EXP, P planes each |
| 2 | IMAGE | LSF\_COUNT | Contribution count, uint16 |
| optional | IMAGE | LSF\_REJ | Rejection counts, low and high, uint16. Off by default (Section 5.4) |
| 3 | BINTABLE | LSF\_SETTINGS | Settings, grid and normalization parameters |
| 4 | BINTABLE | LSF\_WREF | Reference stars of the frame star measures, frozen at initialization, whatever the weighting (Section 10.3) |
| optional | BINTABLE | LSF\_STARS | Cached registration star list (Section 10.4) |
| last | BINTABLE | LSF\_LEDGER | One row per contributing frame |

P is the number of channels: 1 for mono, 3 for OSC (R, G, B).

- **LSF\_ACC** has NAXIS = 3 and NAXIS3 = 4·P. Planes are ordered MEAN, W, M2, W2, each with P channels in R, G, B order. For OSC that is MEAN\_R, MEAN\_G, MEAN\_B, W\_R, ... W2\_B. Each plane is named by an LSFPLn keyword in the LSF\_ACC header, for example LSFPL1 = 'MEAN\_R'; mono files use 'MEAN', 'W', 'M2' and 'W2'. Readers find planes by these names (C6).
  - MEAN is the weighted mean in raw accumulator units, and is the authoritative copy of the mean.
  - W is the sum of weights, M2 the Welford M2 and W2 the sum of squared weights.
  - All four planes are always per channel.
  - **EXP, the exposure map** (optional). It is the sum of EXPTIME, in seconds, over the samples merged at the pixel and channel. When present, NAXIS3 = 5·P, and its planes follow W2, named 'EXP' (mono) or 'EXP\_R', 'EXP\_G' and 'EXP\_B'. NAXIS3 = 5 or 15 never equals the primary's, so stock Siril still sees one image (below). The plane only carries information once exposures differ. While every frame has the same exposure, the exposure per pixel is COUNT times it, so the plane is absent (Section 10.2 says when it is created).
- **Why one cube.** Siril treats a FITS file holding several image HDUs of the same size and BITPIX as a FITS sequence (src/io/fits\_sequence.c, \_find\_hdus). An HDU whose NAXIS3 or BITPIX differs from the first image is skipped. With separate float32 planes, stock Siril opened an LSF file as a 4-frame sequence. NAXIS3 = 4·P or 5·P (4, 5, 12 or 15) never equals the primary's NAXIS3 (none or 3), and LSF\_COUNT and LSF\_REJ differ in BITPIX. Stock Siril therefore sees a single image, without any change to its reader. The same holds for the exposure map, which is why it is a plane of the cube: as a separate float32 HDU of P planes it would match the primary.
- **LSF\_COUNT** has the same grid as the Stack. Mono files use NAXIS = 2. OSC files use NAXIS = 3 with NAXIS3 = 3. LSF\_COUNT may use NAXIS3 = 1 only when every channel always accepts or rejects the same samples (for example cfa.mode = DEBAYER with rejection disabled). It never may for CFADRIZ, where coverage differs per channel.
- **HDU 0 is derived from MEAN.** It is rebuilt from MEAN on every commit: output-normalized when LSFNORM = T, otherwise clamped to \[0, 1\] (Section 4.1), and 0 wherever W = 0. HDU 0 is load-bearing for compatibility (C1), but it is no longer the only copy of the mean, so a damaged primary can be repaired (Section 8).
- **Ledger placement.** An LSF writer puts LSF\_LEDGER last. Readers MUST NOT require it to be last (C6). When rewriting, a writer places unknown HDUs (C7) after the LSF image and settings HDUs and before the ledger, in their original relative order. Optional LSF\_ HDUs, such as LSF\_STARS, also go before the ledger.

**Size.** Per pixel per channel: Stack 4 B, LSF\_ACC 16 B, LSF\_COUNT 2 B, total 22 B. With LSF\_REJ, add 4 B, for 26 B. The exposure map adds 4 B to LSF\_ACC. On a 26-megapixel sensor this is about 575 MB mono and 1.72 GB OSC, or 680 MB and 2.03 GB with rejection maps. Lossless compression of the extensions (Section 5.3) stored them in 23 to 34 percent of their raw size in tests on real OSC stacks. A 7539 x 7302 OSC file of 3.6 GB would shrink to about 1.35 GB, most of which is then the uncompressed primary.

Drizzle at scale s (grid.oversample) multiplies the pixel count by s². A 26-megapixel OSC sensor drizzled at scale 2 gives a 104-megapixel grid: about 6.9 GB uncompressed, of which about 5.6 GB is accumulator and counts. That is more than a writer can usually hold in memory at once (Section 10.10).

## 4. Primary HDU

### 4.1 Image

- BITPIX = -32. NAXIS = 2 (mono) or 3 with NAXIS3 = 3 (OSC, planes R, G, B).
- BSCALE and BZERO MUST be absent or 1 and 0 (Siril writes BSCALE = 1 and BZERO = 0).
- The primary MUST carry ROWORDER, as Siril does ('TOP-DOWN' in the sample header). Every extension image MUST use the same row order, because the registration transforms and any WCS depend on it.
- Pixels with no data (weight 0) are written as 0.0, not NaN. NaN breaks statistics and auto-stretch in many viewers, and 0 is what Siril's output normalization keeps as no data.
- When output normalization is enabled, the primary is MEAN with the same min-max output normalization Siril applies to stacks, and is in \[0, 1\]. Otherwise it is MEAN clamped to \[0, 1\], as Siril clamps its own 32-bit stacks. With the same frames and settings, the primary then equals the image Siril's stack command produces. MEAN is never clamped, so values above 1 (star cores raised by normalization) stay available to the accumulator and to a later output normalization.
- An OSC primary is already RGB data. BAYERPAT, XBAYROFF and YBAYROFF MUST NOT appear in the primary header, because tools that see BAYERPAT may try to debayer an already-debayered image. The pattern is recorded in LSFBAYR instead. This matches Siril 1.4.4, whose OSC stack header carries none of these keywords.

### 4.2 LSF keywords

| Keyword | Type | Meaning |
| --- | --- | --- |
| LSF | logical | T. Marks the file as LSF |
| LSFVER | string | Profile version, e.g. '0.4' |
| LSFUID | string | Unique id, created at initialization, never changed |
| LSFGEN | integer | Generation counter, incremented on every commit, identical in all LSF HDUs |
| LSFSENS | string | 'MONO' or 'CFA' (permanent, Section 1.1) |
| LSFBAYR | string | Bayer pattern (RGGB, GRBG, GBRG, BGGR). CFA only |
| LSFDMOD | string | 'NONE' (mono), 'DEBAYER' or 'CFADRIZ' (permanent, Section 1.1) |
| LSFNORM | logical | T if the Stack is output-normalized |
| LSFNMIN, LSFNMAX | float | Min and max of the output-normalization map |
| LSFHASH | string | SHA-256 of the primary data at last commit, lowercase hex (Section 8) |
| LSFBUSY | logical | Reserved. No writer of this version sets it. A file where it is T is refused for merging (Section 8, check 2) |

LSFHASH is computed over the primary data unit as stored: big-endian float32 values in file order, without the trailing padding to a 2880-byte block. Its 64 characters fit in a string card.

LSFNMIN and LSFNMAX MUST be written with enough digits to round-trip exactly: 9 significant digits for float32 values, 17 for float64. Repair (Section 8) uses them to rebuild the primary from MEAN, and adopt (Section 8) uses them to invert the normalization. A writer MUST refuse to adopt a primary flagged LSFNORM = T if either is missing.

## 5. Extension Image HDUs

### 5.1 Headers

Each extension image carries XTENSION = 'IMAGE', BITPIX, NAXIS and axes, PCOUNT = 0, GCOUNT = 1, EXTNAME, EXTVER = 1, LSFUID and LSFGEN. Cubes (LSF\_ACC, and LSF\_REJ or a per-channel LSF\_COUNT) also carry an LSFPLn keyword naming each plane.

Each also carries **its own copy of the WCS keywords** (CTYPEn, CRVALn, CRPIXn, CD matrix, any distortion terms), but only when the grid has a WCS. INHERIT is not used, because readers support it unevenly. For a cube, the WCS keywords describe axes 1 and 2 only.

DATASUM and CHECKSUM are RECOMMENDED on every HDU. Siril's LSF commands always write them.

### 5.2 Count plane

LSF\_COUNT uses BITPIX = 16 with BSCALE = 1 and BZERO = 32768, the standard unsigned 16-bit convention. The count **saturates at 65535 and never wraps**. A writer that saturates any pixel sets LSFCSAT = T in the LSF\_COUNT header.

Under reliability weights, saturation does not affect n\_eff, which is computed from W2 (Section 10.7). Under noise weighting the count enters the variance estimate as N. There, a saturated pixel's estimate is approximate, with a relative error below 1/65535. Otherwise a saturated count only makes the coverage information at that pixel approximate.

### 5.3 Compression

- **What may be compressed.** Extension images MAY be tile-compressed (ZIMAGE convention). Tables are never compressed.
- **Methods.** Only the methods Siril supports are allowed: RICE\_1, GZIP\_1 and GZIP\_2. HCOMPRESS\_1 and PLIO\_1 are prohibited.
- **Lossless only.**
  - LSF\_ACC uses GZIP\_2 (byte-shuffled, the default) or GZIP\_1, with float quantization disabled (ZQUANTIZ = 'NONE'). RICE\_1 cannot compress floating-point data without quantizing it, so it is prohibited for LSF\_ACC, as is any quantized mode.
  - LSF\_COUNT and LSF\_REJ are integers, which every allowed method compresses losslessly. Siril writes them with RICE\_1, which was nearly as small as GZIP\_2 (14.5 against 13.5 percent of raw) and about 20 times faster in tests.
  - Pre-merge validation refuses any other method, and a quantized LSF\_ACC (Section 8, check 3).
- **Tiling.** The writer uses CFITSIO's default tiles, one row of one plane each. Taller tiles (16 or 64 rows) saved under 1 percent and were slower in tests. One-row tiles also let a writer read and write the accumulator in bands of rows (Section 10.10).
- **Primary.** The primary HDU MUST NOT be tile-compressed, and MUST NOT be replaced by an empty primary followed by a compressed image extension. That pattern breaks C1. Siril's existing compressed save writes exactly this pattern, so the LSF writer cannot use it for HDU 0. It enables compression only after HDU 0 is written, and disables it before the tables.
- **Not stored in the file.** Whether to compress is a choice made when writing, not a fact a reader needs, so LSF\_SETTINGS does not store it. Whether the file is compressed, and how, is read from the extensions themselves (ZCMPTYPE of LSF\_ACC).
- **Commands.**
  - lsf\_init writes uncompressed extensions by default. -compress uses GZIP\_2 for LSF\_ACC, -compress=gzip1 uses GZIP\_1.
  - lsf\_append keeps the compression of the file, unless -compress or -nocompress is given.
  - lsf\_util -compress and -decompress convert a file without appending, replacing it through the commit protocol or, with -out=, into a new file (Section 8).
  - Siril's global FITS compression preference does not apply to LSF files, because it allows quantized float compression. It does apply to the plain images that lsf\_export and lsf\_maps save (Section 10.8).
- **Cost.** On a 7539 x 7302 grid, one channel's LSF\_ACC and LSF\_COUNT took about 4.5 s to write compressed instead of 1.6 s, and 2.0 s to read instead of 1.0 s (from the operating system's cache). An OSC commit therefore takes about 9 s longer. CFITSIO compresses on one thread, and its gzip level cannot be set through its public interface (CFITSIO 4.7).
- **Whole-file gzip** (.fits.gz) is allowed for transport only. A file must be decompressed before appending.
- **Updates.** Every append rewrites the whole file (Section 10.2), so compressed HDUs never need to be updated in place.

### 5.4 Rejection planes (optional)

LSF\_REJ is off by default. It is created only when the user asks for it, by an option on lsf\_init or lsf\_append. It records, per pixel and channel, how many samples were rejected low and how many high.

- It uses BITPIX = 16 with BSCALE = 1 and BZERO = 32768.
- NAXIS3 = 2·P, ordered REJLO (R, G, B) then REJHI (R, G, B), with planes named by LSFPLn ('REJLO\_R' and so on; 'REJLO' and 'REJHI' for mono).
- Low and high are always stored separately, never merged. The two directions mean different things: high rejections come from satellites, aircraft, hot pixels and cosmic rays, low rejections from cloud, cold pixels and changing dust shadows.
- It counts rejections from both the batch-local test and the history test (Section 10.6).
- Counts saturate at 65535 and never wrap. A writer that saturates any pixel sets LSFRSAT = T in the LSF\_REJ header.
- When rejection maps are enabled on a later append rather than at initialization, they count only from that commit on. LSFRGEN in the LSF\_REJ header records the first generation counted.
- The setting rejmaps.enabled mirrors whether the HDU exists.

## 6. LSF\_SETTINGS

A BINTABLE with one row per setting: KEY (32A), TYPE (8A: STR, INT, FLT, BOOL or JSON) and VALUE (1PA, variable-length string). A table is used instead of header keywords because settings include nested structures (rejection parameters, alignment options), and new settings should not need a schema change.

Required keys:

- sensor.class, sensor.bayer, sensor.binning, cfa.mode
- grid.oversample (the drizzle scale, 1 without drizzle)
  - The grid size and WCS are the primary's NAXIS and WCS keywords, and are not duplicated.
- grid.origin\_x, grid.origin\_y (margins of the grid before the first reference frame, from max framing at initialization and from grid expansions, Sections 10.1 and 10.5)
- grid.undistorted (BOOL): every frame was undistorted before it was merged, so the grid is a plain projection of the sky (Section 10.4). An absent row means F
- align.transform, align.interpolation, align.min\_pairs, align.findstar (JSON: the star-detection parameters used for registration)
- align.disto (JSON: source, 'file', 'master' or 'frames', and file, the model file or the master pattern, empty for 'frames'), the distortion model the appends undistort their frames with. 'frames' is the SIP terms of each frame's own WCS, as Siril's astrometric registration undistorts them. Required when grid.undistorted is T, absent otherwise
- reject.method, reject.params, reject.hist\_low, reject.hist\_high, reject.hist\_margin (Section 10.6)
- weight.scheme, plus the frozen weight constants for that scheme (Section 10.3):
  - noise: weight.noise\_norm
  - wfwhm: weight.fwhm\_min, weight.fwhm\_max, weight.fwhm\_norm
  - nbstars: weight.stars\_min, weight.stars\_max, weight.stars\_norm
  - wfwhm and nbstars: weight.model ('range' or 'inverse'), weight.floor and weight.range (FLT). Files written before them have neither row, and mean range, 0.01 and 1
  - A drizzled file, whatever the scheme: weight.driz\_norm (JSON, one value per channel). CFADRIZ says so for a CFA sensor; on a mono sensor, where cfa.mode is NONE, this key does
- norm.method, norm.estimator ('ikss' or 'lite', the statistics the normalization uses), norm.rgb\_equal (BOOL)
  - norm.baseline (JSON: loc and scale, one value per channel each), the reference statistics frozen from the first batch. Required unless norm.method is none
- outnorm.enabled, outnorm.min, outnorm.max
- rejmaps.enabled
- filter (mono and OSC files), created\_utc, writer (name and version)

Conditional key: ledger.orig\_baseline, the baseline header that ORIGHDR diffs are taken against (Section 9). Required when any ledger row stores a diff. A writer MAY store diffs. Siril stores every header verbatim and writes no baseline.

**Facts and defaults.** The rows are of two kinds.

- **Facts** are values the data depends on, which a reader or a later append needs to interpret the accumulator: sensor.\*, cfa.mode, grid.\* (grid.undistorted included), weight.\*, norm.\* and outnorm.\*. They never change after initialization, except outnorm.min and outnorm.max, which change on every commit, grid.origin\_\*, which changes on grid expansion, and weight.model with weight.fwhm\_norm or weight.stars\_norm, which change once when an append switches the weight model to inverse (Section 10.3).
- **Defaults** are the processing choices of the first batch (align.\*, reject.\*). They are recalled on append, so later batches are processed the same way unless the user overrides them for that append. An overridden align.disto replaces the stored one, as an overridden rejection does.

Choices that apply to one append only are command options and are not stored. These are whether and by what step the grid expands (Section 10.5), the registration method (Section 10.4), normalization on the overlap with the Stack and feathering (Section 10.11), and the drizzle pixel fraction and kernel. The ledger records what each frame used. An append to a CFADRIZ file uses the pixel fraction and kernel of the last ledger row unless told otherwise, and another one needs lsf\_append -drizzle\_force (Section 10.3). Files from earlier writers may hold grid.reference, grid.pixfrac, grid.kernel, grid.expand and grid.expand\_step rows; a reader ignores them and a writer drops them.

Defaults are recalled on append. Some values are mirrored in the primary header: sensor.class and LSFSENS, sensor.bayer and LSFBAYR, cfa.mode and LSFDMOD, outnorm.\* and LSFNORM, LSFNMIN, LSFNMAX, and rejmaps.enabled and the presence of LSF\_REJ. These MUST agree. A disagreement fails validation (Section 8).

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
| HISTORY | One card per append, giving the generation, batch number and frame count, followed by one card per channel with the percentages of samples rejected low and high by the batch-local and the history tests. Repair, Adopt (Section 8), each change of the compression of the extensions (Section 5.3), each grid expansion (Section 10.5) and each WCS import (Section 8) add one card each |
| WCS keywords | Updated on grid expansion (CRPIX, Section 10.5), and replaced when a plate-solved WCS is imported (lsf\_util -wcs, Section 8) |
| MIPS-HI, MIPS-LO | Removed on append, because display cutoffs saved for the previous stack no longer fit |

**Set at initialization, then checked on append.** A mismatch fails the compatibility check and the frame is not merged:

- XBINNING, YBINNING
- XPIXSZ, YPIXSZ and FOCALLEN, scaled for drizzle the way Siril does it. Raw frames carry the sensor's own values, so an incoming frame's values are scaled the same way before they are compared. FOCALLEN matches within 5 percent: Siril's plate solving replaces it with the focal length it measures, a few percent from the nominal one, while a reducer or another telescope differs by much more
- INSTRUME, FILTER

**Copied from the first frame.** If a later frame differs, the keyword is removed from the primary, and the ledger keeps the per-frame values:

- EXPTIME (the sub-exposure length). LIVETIME stays exact when it is removed.
- GAIN, OFFSET, CCD-TEMP, EGAIN
- TELESCOP, OBJECT

**Never in the primary:** BAYERPAT, XBAYROFF, YBAYROFF (Section 4.1).

## 8. Integrity Binding

The primary is a view of MEAN. A foreign tool that edits it while preserving the extensions would produce a file whose image no longer matches its accumulators. FITS has no way to prevent this. The profile makes it **detectable** and **repairable**.

**Binding:** LSFUID and LSFGEN in every LSF HDU, LSFHASH in the primary, and DATASUM and CHECKSUM on every HDU.

LSFHASH binds the primary instead of DATASUM, because DATASUM is a ones-complement sum of 32-bit words, and its result does not depend on word order. A tool that flips, mirrors or rotates the primary by 180 degrees leaves DATASUM unchanged, and a DATASUM comparison would accept the file. SHA-256 depends on position. DATASUM and CHECKSUM stay as the standard transport checks.

**What the hash does.** It detects; it cannot prevent. It cannot stop another program overwriting the file. It makes Siril's LSF commands refuse to merge into a primary that no longer matches MEAN. Within Siril, the LSF commands are the only path that writes an LSF file, and every write goes through the commit protocol (Section 10.2).

**Pre-merge validation.** A writer MUST check, in order:

1. LSF = T and LSFVER is supported (Section 11).
2. LSFBUSY is absent or F.
3. Every required HDU is present (PRIMARY, LSF\_ACC, LSF\_COUNT, LSF\_SETTINGS, LSF\_WREF, LSF\_LEDGER, found by EXTNAME) with consistent dimensions, LSF\_ACC holds the expected named planes (the EXP planes too when NAXIS3 = 5·P), and any tile compression of LSF\_ACC and LSF\_COUNT uses an allowed method losslessly (Section 5.3). LSF\_REJ, when present, is checked the same way. LSF\_WREF holds at least one star.
4. LSFUID and LSFGEN are identical across all LSF HDUs (LSF\_REJ included when present).
5. DATASUM and CHECKSUM match on every extension where they are present.
6. Every ledger GEN is at most LSFGEN. A generation with no ledger rows is one that merged no frame (Repair, Adopt, a change of compression or a WCS import), and is recorded in a HISTORY card.
7. Mirrored settings agree.
8. The SHA-256 recomputed from the primary data equals LSFHASH.

If any of checks 1 to 7 fails, the writer refuses to merge, reports which check failed, and leaves the file untouched. The file stays readable as an image (C1). There is no silent repair. If LSF\_ACC itself is damaged (check 3 or 5), the accumulator cannot be trusted, and the user restores a backup copy.

If only check 8 fails, the accumulators are intact but the primary was changed outside Siril. The writer stops, reports that the Stack was modified outside Siril, and offers two explicit actions. In Siril they are options of lsf\_util, whose -check option only validates the file and summarizes it:

- **Repair** (recommended, lsf\_util -repair). Rebuild HDU 0 from MEAN (Section 3). This undoes the outside edit: the rebuilt primary is bit-identical to the one written at the last commit.
- **Adopt** (lsf\_util -adopt). Accept the edited primary as the new mean. A pixel whose primary value is still the one its MEAN gives was not edited and keeps MEAN exactly, which also keeps values the primary clamps (Section 4.1). An edited pixel gets MEAN from HDU 0, inverting output normalization with LSFNMIN and LSFNMAX when LSFNORM = T. Pixels with W = 0 stay 0. M2 is not adjusted, so adopt is meant for small cosmetic fixes the user wants to keep. Adopt cannot work if the dimensions have changed. The primary is then rebuilt from the new MEAN, so with output normalization an edit beyond the old range rescales the whole Stack.

Either action records a HISTORY card, increments LSFGEN, and recomputes LSFHASH and the checksums.

lsf\_util has two more actions:

- **Compression.** -compress\[=gzip2|gzip1\] or -decompress changes the compression of the extensions (Section 5.3).
- **WCS import.** -wcs=file imports a WCS. The file is an image of the Stack that was plate-solved, such as the output of lsf\_export (Section 10.8) after plate solving, or a Siril .wcs file (NAXIS = 0, with IMAGEW and IMAGEH for the dimensions). It must have the Stack's dimensions. When it carries the export's provenance HISTORY card, that card must name the file's LSFUID. The writer replaces the WCS keywords of HDU 0 and of every extension (Section 5.1). The pixels do not change, so neither does LSFHASH.
  - Siril's plate solving may flip an image whose solution shows it upside down, and it then swaps its ROWORDER. When the image's ROWORDER differs from the Stack's, its WCS describes the flipped pixels, and the writer flips it back.
  - On an undistorted grid (grid.undistorted = T), the writer drops the SIP terms of the imported WCS and logs the largest correction they made at the edges of the grid, with a warning above 0.5 px. The grid is a plain projection by construction, so these terms only fit residuals, and Siril would distort the Stack's stars with them when it registers against it.
  - A WCS the file already has, to 1e-4 arcsec at the corners, the edges and the centre, is reported and changes nothing.

Both actions record a HISTORY card and increment LSFGEN although no data changes, so that two different files never share an LSFUID and LSFGEN. lsf\_util needs at least one action. Its actions combine into one rewrite with one LSFGEN increment, for example -repair -compress. Every rewrite rebuilds HDU 0 from MEAN. A change of compression or a WCS import without -repair or -adopt is therefore refused while check 8 fails, because it would silently repair the Stack. A requested state the file already has is reported and changes nothing. With -out=, the result is written to a new file and the source is left untouched; with no other action this makes a copy. The copy keeps the LSFUID of its source with an incremented LSFGEN and a HISTORY card naming the source, so the two files are successive versions of one stack, and only one of them should be appended to.

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

These are design expectations, not observed behavior, except for the Siril behavior in the second row, which was confirmed against the Siril source.

## 9. LSF\_LEDGER

A BINTABLE, written as the last HDU (Section 3), with one row per contributing frame, carrying its identity, acquisition data, registration, normalization, weight and quality measures. Per-plane columns have repeat count P. Their TFORM is written below as nD or nE, where n = P: 1D for mono, 3D for OSC.

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
| ALIGN | 1PA | Alignment parameters used (JSON). On append, method is 'stars' or 'wcs' (Section 10.4); disto is 'file', 'master' or 'frames' when the frame was undistorted |
| REG\_DIST | 1PA | The distortion model the frame was undistorted with: CRPIX1, CRPIX2 and the SIP cards (A\_ORDER, A\_p\_q, B\_..., AP\_..., BP\_...), stored as 80-character cards like ORIGHDR. Empty when the frame was not undistorted |
| DRZ\_PIXF | E | Drizzle pixel fraction. Only in drizzled files |
| DRZ\_KERN | 8A | Drizzle kernel: square, point, turbo, gaussian, lanczos2 or lanczos3. Only in drizzled files |
| NORM\_SC, NORM\_OFF | nD | Normalization scale and offset per plane |
| NORM\_MTH | (12·P)A, TDIM = (12,P) | linear, offset\_only, failed or skipped, per plane (12A mono, 36A OSC). Siril writes norm.method, with '+ovl' when the frame was normalized on the overlap with the Stack (Section 10.11) |
| FEATHER | J | Feathering distance the frame was merged with, px, 0 for none (Section 10.11). Added by the first commit that feathers; older rows get 0 |
| WEIGHT | nD | Weight applied per plane, after the scaling of Section 10.3 |
| FWHM, WFWHM, ROUND, REGQ, REGBG | E | Quality metrics |
| NSTARS | J | Star count |
| BKG\_MED, NOISE, NOISE\_N | nE | Per-plane background and noise |
| STAT\_MEAN, STAT\_SIG, STAT\_MAD, STAT\_SQBW, STAT\_LOC, STAT\_SCL | nE | Per-plane Siril image statistics |

Rules:

- Readers locate columns by TTYPE. Every quality column is optional. Unknown values are NaN (floats), TNULL (integers) or empty (strings). Files written by earlier versions simply lack columns.
- Unknown columns MUST be preserved on rewrite (C7).
- String columns hold printable ASCII only (32 to 126), as FITS requires. ORIGHDR is stored the way FITS stores a header: 80-character cards one after another, padded with spaces, with no line breaks. A reader splits it every 80 characters. Siril keeps headers as text with one card per line, so the writer removes the line breaks.
- The ledger is append-only in normal operation. Because the writer puts it last, the writer can rewrite it in full on every append without moving any other HDU. It holds scalars and short text, so this stays cheap.
- REG\_XFRM always maps onto the grid as it is now. A grid expansion rewrites the existing rows (Section 10.5). It is exact, because every registration this profile supports is a homography (Section 10.4). When REG\_DIST is set, REG\_XFRM maps the frame's undistorted pixels: the undistortion applies on the frame's side and the expansion shift on the grid's, so the shift stays exact.

## 10. Operations

### 10.1 Initialize

lsf\_init does the following:

1. Builds the canonical grid from the first batch, using Siril's existing registration. When the registration is not applied yet, the writer applies it. By default the grid is the union of all the first batch's frames (max framing), so no frame is cropped, as grid expansion does on append (Section 10.5). The user may instead ask for Siril's other framings: the reference frame's own size (current), the area common to all the frames (min), or the reference frame's size centred on the frames' centre of gravity (cog). With max framing every frame is placed on the whole grid, rather than cropped to its own bounding box as Siril's max framing does. grid.origin\_x and grid.origin\_y then hold the margins of the grid before the reference frame (negative when the grid starts inside it, as with min framing), and REG\_XFRM maps each frame onto the grid. A batch whose registration is already applied defines the grid by the size of its frames. A batch registered with distortion correction (Siril's register -disto=image, file or master) is undistorted by Siril's apply, which also drops the SIP terms from the WCS it gives the frames: the grid is then undistorted for good (grid.undistorted = T, align.disto the model, with a model file stored as an absolute path; Section 10.4). A batch registered through plate solving (Siril's astrometric registration, seqplatesolve without -noreg) is applied the same way, and can hold several mosaic panels (Section 10.11). The writer computes the transforms from the WCS of the frames with Siril's functions, for the framing asked. Max framing uses Siril's projection at the centre of the frames, at the angle that frames them best, whichever way the reference frame is turned; the grid is then the union of the frames on it, and grid.origin\_\* the margins before the reference frame's corners, however it lies. Siril's apply keeps these transforms and undistorts each frame with the SIP terms of its own WCS. When some frame has SIP terms, the grid is undistorted (grid.undistorted = T, align.disto 'frames'); when none has, it is not. A sequence Siril has drizzled is refused, on append as well: its weights may be rescaled and reduced to 8 bits, and its drizzle parameters are not recorded. CFA frames that are not debayered are refused too, since they would be merged as a mono image of the raw mosaic, unless the user asks for a drizzled file (lsf\_init -drizzle, with -scale=, -pixfrac= and -kernel=). The batch is then a calibrated CFA sequence that is not debayered, or a mono one, with a registration that is not applied yet. The writer applies its registration with Siril's drizzle at scale grid.oversample, a whole number from 1 to 3, which gives frames and their drizzle weights (Section 10.3), RGB for a CFA sequence. The framing is chosen in the frames' pixels, so the grid and grid.origin\_x and grid.origin\_y are multiples of the scale. REG\_XFRM maps onto the grid itself. Siril's drizzle multiplies the values by the square of the scale, so the Stack, which is clamped to \[0, 1\] (Section 4.1), clips far sooner unless it is output-normalized; MEAN keeps them. This matches Siril's stack. Feathering, -overlap\_norm and a registration with distortion correction or through plate solving are not supported with it yet, nor are lsf\_append -noreg and -wcsreg on such a file.
2. Merges the first batch into the empty accumulator, as an append does (Section 10.2, steps 2 to 5).
3. Freezes the weight constants, including weight.driz\_norm under CFADRIZ (Section 10.3).
4. Writes every HDU. LSF\_REJ is written only if the user asked for rejection maps.
5. Sets a new LSFUID and LSFGEN = 1, and writes LSFHASH and the checksums.

If the first batch has fewer than 25 frames, the writer SHOULD warn the user before creating the file. Per-pixel scatter estimated from so few frames is unreliable, and any rejection decision made on it stays in the stack. The warning does not block creation. The profile defines no minimum stack size for rejection, so rejection behaves the same at every count. The threshold of 25 is a default.

### 10.2 Append (commit protocol)

1. Run the Section 8 validation.
2. Process the batch with Siril's existing code:
   - the compatibility check: same channel count, sensor class and Bayer pattern, and the Section 7 keywords set at initialization (binning, pixel size, focal length, instrument, filter). A frame that differs is not merged;
   - calibration;
   - registration onto the grid (Section 10.4);
   - normalization with Siril's formulas against the frozen baseline (Section 1.1).
3. Grow the grid first if the registered frames reach beyond it, unless the user asked to crop them to the grid (Section 10.5).
4. Compute the frame weights (Section 10.3) and run rejection (Section 10.6).
5. Update the accumulator, in memory or in temporary files (Section 10.10).
   - The batch is merged into history with the weighted pairwise update. With δ = MEAN\_b − MEAN\_a and W = W\_a + W\_b:
     - MEAN = MEAN\_a + δ·W\_b / W
     - M2 = M2\_a + M2\_b + δ²·W\_a·W\_b / W
     - W2 = W2\_a + W2\_b
     - COUNT = COUNT\_a + COUNT\_b, saturating
     - EXP = EXP\_a + the sum of EXPTIME over the batch's merged samples, when the file has an exposure map
   - LSF\_REJ is updated if present.
   - **The exposure map is created** by the first commit that merges a frame whose known exposure differs from the others' (at initialization, when the first batch is already mixed). The writer does not need to ask for it, and logs it with a HISTORY card. On an append the existing pixels get EXP = COUNT × the exposure every merged frame had, which is exact, and only a lower bound where a count saturated (LSFCSAT). That needs every ledger row to have a known exposure, all the same. Otherwise (a frame of unknown exposure, or a file that mixed exposures before exposure maps existed) the writer warns that the file cannot have an exposure map, and continues without it.
   - **A frame of unknown exposure** in a file that has an exposure map is merged, and adds 0 s to EXP. The writer warns, naming the frames, and a HISTORY card gives their number. The map is then a lower bound wherever they contributed.
6. Rebuild HDU 0 from MEAN, update the Section 7 keywords and append the ledger rows. Increment LSFGEN in every LSF HDU, and recompute LSFHASH and the checksums.
7. Write the complete new file to a temporary name in the same directory.
8. Flush to disk, then atomically replace the original (Section 10.10). A writer SHOULD offer to write the result to a new file instead, leaving the original untouched (lsf\_append -out=), or to keep the previous file as a backup (-backup). As with lsf\_util -out= (Section 8), the new file keeps the LSFUID of its source, with the incremented LSFGEN and a HISTORY card naming the source, and only one of the two files should be appended to.

### 10.3 Frame weights

Siril's existing weighting functions normalize the weights of one stacking run so that their mean is 1. The wFWHM and star-count weights are also scaled between that run's minimum and maximum, so the worst frame of a run gets weight 0. If these were applied to each batch independently, every batch would carry the same weight per frame whatever its quality, and weights would not be comparable across appends.

LSF therefore reuses Siril's per-frame inputs and Siril's formulas, but with constants frozen from the first batch:

- **noise:** the raw weight is 1 / (pscale² · bgnoise²), using Siril's normalization scale and background noise for the channel. weight.noise\_norm is the mean raw weight of the first batch, and every later weight is the raw weight divided by it. No clamping is needed.
- **wfwhm:** Siril's formula, with weight.fwhm\_min and weight.fwhm\_max taken from the first batch, and weight.fwhm\_norm the first batch's mean weight. A later frame's wFWHM is clamped to \[fwhm\_min, fwhm\_max\] before the formula is applied, and the result is raised to the weight floor (below).
- **nbstars:** Siril's formula, with weight.stars\_min, weight.stars\_max and weight.stars\_norm from the first batch. Star counts are clamped the same way, and the result is raised to the weight floor.
- **Range span** (weight.range, lsf\_init -weight\_range=, 1 by default). The span of the first batch is multiplied by it before it is frozen: weight.fwhm\_max = fwhm\_min + range · (worst − fwhm\_min), and weight.stars\_min = stars\_max − range · (stars\_max − fewest), at least 0. 1 is Siril's rule. A larger value leaves room for later frames worse than the first batch's worst.
- **Weight floor** (weight.floor, lsf\_init -weight\_floor=, 0.01 by default, between 0 and 1). Siril's wfwhm and nbstars formulas give 0 to the worst frame of a run. With constants frozen from the first batch, every later frame at or beyond that frame's quality would get 0 and be dropped, which can remove a whole night of slightly worse data. The formula result, before division by the norm, is therefore at least the floor, where 1 is the best frame of the first batch. It applies to the first batch too, so its norm includes it.
- **Weight model** (weight.model, for wfwhm and nbstars). Every file starts with **range**, the formulas above, so that the first Stack follows Siril's rules. An append can switch it to **inverse** (lsf\_append -weight\_model=inverse), which has no range: the raw weight is (weight.fwhm\_min / wFWHM)², Siril's formula with an infinite worst, or (stars / weight.stars\_max)², raised to the floor. Nothing is clamped any more, and the first batch's constants become a common scale of every weight.
  - The switch applies from that append on. The frames already merged keep their weights, which the accumulator holds and cannot change frame by frame. The new norm is chosen from the ledger so that, under inverse, those frames would weigh in total what they weigh now: the new frames are on the same scale on average, and the accumulator is not rescaled.
  - The switch is one way: going back to a range would clamp frames that inverse weighted freely.
- **The range check on append.** After the star measures of a batch, the writer counts its frames beyond the frozen range: better than the first batch's best, which all get the top weight whatever their differences, and worse than its worst, which all get the worst's weight or the floor. When there are any, it warns and recommends the switch to inverse. A frame's weight then no longer follows its quality.

**Why weighting a living stack is harder than weighting one run.** These rules follow from what an accumulator can and cannot do:

- **Weights are relative.** A weighted mean only uses the ratio of one frame's weight to the others'. Siril sees all the frames of a run before it weights any of them; a Living Stack must weight each batch as it arrives, against frames it will never see again. The constants that set the scale (the best and worst of the first batch, or the mean noise) are therefore frozen at initialization.
- **The past cannot be re-weighted, except by one common factor.** The accumulator keeps only sums (MEAN, W, M2, W2), not the frames. Multiplying every past weight by the same factor is exact: W becomes c·W, W2 c²·W2, M2 c·M2, and MEAN does not change. Changing the ratios between past frames, as a new range would, needs the frames themselves. This is why a change of range, model or scheme can only apply going forward. A file that should be weighted differently from the start must be made again from the frames with lsf\_init.
- **A range frozen from a batch that is not typical clamps later frames.** Under range, a first batch of the best frames sends later, worse ones to the floor, and a first batch of the worst gives every later, better frame the same top weight. On 499 frames, best first lost 20 percent SNR, with 37 percent of the frames at the floor; random first batches stayed within a few percent of Siril, with about 1 percent at the floor. Inverse weights do not depend on the first batch, at the cost of weaker contrast: the worst frames keep about a quarter of the best one's weight, where range takes them to the floor. That stack was about half a percent softer than Siril's, with 1 to 4 percent more SNR.
- **Only wfwhm and nbstars have a range.** noise weights are already the raw weight over a frozen mean, a common scale like inverse, and none has nothing to scale. The range options, the check and the switch therefore apply to wfwhm and nbstars only.
- **A metric needs the same reference every time.** wFWHM counts the reference stars a frame misses, so the reference must be the same for every frame of every batch, or later frames are measured against a different, usually deeper, list (see LSF\_WREF below).
- **The ledger records what each frame got, so what any other scheme would have given can be computed.** Applying it needs the frames.
- **none:** every weight is 1.

**Drizzle pixel weights.** In a drizzled file, Siril's drizzle gives each output pixel of a frame a weight per channel, 0 where no input pixel of that channel's color landed. Siril divides these weights per channel by max\_weight, the largest weight of a small uniform patch drizzled with no transform at the scale, pixel fraction, kernel and CFA pattern of the run, and by default saves them as 8-bit integers, which saturate wherever a real frame's weight exceeds max\_weight (at scale 2, with any pixel fraction below 1). LSF therefore does the following:

- It keeps the drizzle weights in float32 and unscaled.
- It divides them by weight.driz\_norm, the per-channel max\_weight of the first batch, frozen at initialization. The result is d, which can exceed 1.
- A sample's weight is the frame weight of its channel times d.
- A sample with d = 0 is absent: it is not counted, tested or rejected. So is one with a negative d, which the lanczos kernels can give, where Siril's stack keeps it.
- W, W2, COUNT, LSF\_REJ and the history test (Section 10.6) all use this per-sample weight.

With the frozen norm, the total weight of a frame does not depend on its pixel fraction under the square kernel: each input pixel's weights sum to 1, so a smaller drop puts more weight on fewer pixels. Between kernels it holds at scale 1 only: at scale 2 a lanczos frame carries 4 times the weight of a square one, and a turbo frame 16 times.

The frames of a file should keep one pixel fraction and kernel, and lsf\_append refuses others without -drizzle\_force. A frame's noise moves with its drop: at a pixel fraction of 0.5 against 1, on real and synthetic frames, the background noise rises by 15 to 25 percent, and the normalization scale of a mono frame by 12 to 17 percent, of a CFA frame's green by 2 to 13 percent, with red and blue unchanged since their drops rarely overlap. Under a scaling normalization that would scale the signal of those frames down, and shift the colour of a CFA frame.

**Feathering weights.** With feathering (Section 10.11), a sample's merge weight is also multiplied by its ramp m, between 0 and 1 near the frame's edges. A sample with m = 0 is absent. The history test predicts the spread with the frame weight alone, because the ramp is a blending choice and says nothing about the sample's noise. Under noise weighting, σ0² = M2 / (N − 1) assumes weights proportional to 1/σ², so it is biased low inside the feather band. Reliability weights are not affected.

**Star measures against frozen reference stars (LSF\_WREF).** Siril measures a frame's wFWHM when it registers it: the FWHM, raised by the share of the reference frame's stars that found no pair in the frame (2·FWHM·(fitted − pairs)/fitted + FWHM). Within one Siril run the reference is one frame. Across appends it would not be: a first batch registered within itself pairs nearly all its stars, while later batches registered against the Stack, which has more stars than any frame, cannot. Their wFWHM then grows with the depth of the Stack, not with their quality, and with frozen constants every appended frame ends at the weight floor (seen in tests: every frame of every append). So LSF measures every frame itself, the same way, whatever registered it:

- At initialization, whatever the weighting, the stars of the reference frame of the first batch, on the grid, become LSF\_WREF, a BINTABLE with columns X and Y in FITS pixels of the grid (1-based), carrying LSFUID and LSFGEN. They never change, except that a grid expansion shifts them as it shifts CRPIX (Section 10.5).
- Each frame to merge is measured on the grid, after its registration and any undistortion: its stars are detected with align.findstar (in a drizzled file, on a copy whose holes, the pixels no input pixel reached, are interpolated along their row and column, and which at a scale above 1 is first binned back to the frames' scale and level, the mean of the data in each block, since a frame's pixels and holes are that wide on the grid; positions and FWHM are then given in grid pixels), its FWHM and roundness come from them, and a pair is a frame star within 3 pixels of the frames (the inlier radius of Siril's star alignment) of a reference star, each reference star used once. fitted is the number of reference stars where the frame has data, capped at 2000 as star alignment caps it, so that a frame framed differently is not counted as worse.
- Its wFWHM is 2·FWHM·(fitted − pairs)/fitted + FWHM, and its star count the number of stars detected. The wfwhm and nbstars weights use these, and the ledger's FWHM, WFWHM, ROUND and NSTARS columns record them.
- The first batch is measured the same way, so its frozen constants (fwhm\_min, fwhm\_max and the norm) are on the same scale as every later frame. The reference frame's own wFWHM is its FWHM, as in Siril.
- LSF\_WREF is required. Under noise or no weighting the weights do not use it, but the ledger's star columns stay comparable across batches, for frame selection and for any later change of metric. Its cost is one star detection per frame; the table itself is small (at most 2000 stars, about 37 KB). A reference frame with no star refuses the initialization.

On the same Stack, frames registered with stars or through their WCS got the same pair counts, and repeated appends of the same frames the same wFWHM to about 1 percent.

A frame whose weight is still 0 contributes nothing: one without a wFWHM (no stars measured), or with a zero noise weight. It is not merged, and it is reported to the user. These rules live in an LSF helper alongside Siril's code; Siril's own weighting functions are not changed. Siril's nbstack weighting is not supported in this version.

### 10.4 Registration and the grid reference

The canonical grid is the Stack's own pixel grid. At initialization it is defined by the first batch, by default as the union of its registered frames (Section 10.1), scaled by grid.oversample. That frame is not stored, so later batches are registered against **HDU 0 itself**, by one of two methods chosen for each append. The choice is a command option and is not stored (Section 6).

- **Stars** (the default). Siril's existing star detection, with the parameters stored in align.findstar, and Siril's existing global registration against HDU 0. The Stack is the highest-SNR reference available and already defines the grid.
- **WCS** (lsf\_append -wcsreg). Each frame's transform onto the grid is computed from the frame's WCS and the Stack's WCS, as Siril's astrometric registration does, but fitted with the transform type of align.transform on an 8 by 8 grid of points over the frame instead of its four corners. The frames must be plate-solved; a frame that is not is not merged. The Stack must have a WCS (Section 8, WCS import). This method suits frames that share too few stars with the Stack. The frames' star measures, for the quality columns and the wfwhm and nbstars weights, are made on the grid against LSF\_WREF, as for every other frame (Section 10.3), so frames are weighted alike whichever method registered them. A registration made beforehand against the file (register -extref=) is not used.

Every transform is a homography, so REG\_XFRM describes it exactly (Section 9).

**Distortion correction** is a property of the grid, chosen at initialization (Section 10.1) and recorded in grid.undistorted:

- **An undistorted grid** undistorts every batch before merging it, with the model of align.disto or the one given to lsf\_append with -disto= (as Siril's register takes it: master, or file and its name; or frames), which then replaces it. The frames' stars are undistorted with the model before they are matched, and Siril's apply resamples the frames through the model and the homography. The Stack itself is never undistorted as a reference, since it is undistorted already. Under -wcsreg, the transform maps the linear part of the frame's WCS onto the Stack's linear WCS: between two TAN projections this map is exactly a homography. REG\_DIST records the SIP terms each frame was undistorted with.
- **A grid without distortion correction** never undistorts. Its Stack keeps the distortion of its first batch, and frames undistorted later would not line up with it toward the edges, so -disto= is refused, as is a registration against it made with distortion correction.
- An existing file cannot switch. lsf\_append -noreg on an undistorted grid merges the frames as they are, with a warning: they must be undistorted already.

- **The frames' own SIP terms** (align.disto 'frames', from an initialization by plate solving, or lsf\_append -disto=frames). Siril's star alignment undistorts the stars with one model, so these frames are registered only through the WCS (-wcsreg), and must be plate-solved; star registration is refused unless one model is given with -disto=file or master. Siril's apply then undistorts each frame with the SIP terms of its WCS while keeping the transforms onto the Stack: for an external reference it now takes only the distortion from the frames' astrometry, without computing the transforms again against a frame of the sequence. REG\_DIST records each frame's own terms. A frame without SIP terms is not undistorted.

In a drizzled file the incoming frames are at the sensor's scale, and CFA images under CFADRIZ, while the Stack is at scale grid.oversample, and RGB under CFADRIZ. The transforms are computed onto the grid at its scale, so the drizzle step must not apply the scale a second time: the apply gets S⁻¹·H, with S = \[s 0 t; 0 s t; 0 0 1\] and t = (s − 1) / 2, the scaling about the pixel centres that it applies itself. Plain division by s leaves a shift of t. REG\_XFRM keeps H. Stars are detected on the Stack as they are, which in tests at scales 2 and 3 placed the frames within 0.15 grid pixels rms of a registration on the frames' own scale.

The optional LSF\_STARS BINTABLE caches the star list detected on HDU 0, typically under 1 MB, to make appends faster and repeatable. Its columns are X, Y, FLUX and FWHM, in grid pixels. Its header carries LSFGEN. A reader uses the cache only when that LSFGEN equals the file's; otherwise it detects stars again.

### 10.5 Grid expansion

Framing drifts between sessions, so the grid must be able to grow. When registered frames would cover pixels outside the current grid, the writer grows the grid before merging, unless the user asked to crop the batch to the grid instead (current framing). Expanding is the default, as max framing is at initialization (Section 10.1). This choice and the step below belong to the append command or the user interface, not to the file: they change no data a reader interprets, so LSF\_SETTINGS does not store them. On append only max and current framing apply: min and cog would shrink or move the grid. Older files may hold grid.expand and grid.expand\_step rows, which are ignored (Section 6).

A pixel outside the grid counts as covered when a transformed frame reaches its center, so sub-pixel overhangs from registration noise need no expansion. Registration transforms (REG\_XFRM) use Siril's registration coordinates, whose y axis runs opposite to the FITS storage order: rows a frame needs below y = 0 there are inserted after the existing data, and rows it needs past the last row are inserted before it. Then:

- Each side grows in multiples of the expansion step (default 64 pixels), so small drifts do not rewrite the file on every append.
- Every image plane is padded with zeros, the exposure map included. W = 0 means no data, so padding changes no existing pixel.
- Rows or columns inserted before the existing data, in storage order, are added to CRPIX1 and CRPIX2 in every HDU that carries a WCS, and to grid.origin\_x and grid.origin\_y.
- REG\_XFRM in every existing ledger row is rewritten to include the shift, so every transform keeps mapping onto the current grid.
- LSF\_STARS is discarded.
- The X and Y of LSF\_WREF get the same shifts as CRPIX1 and CRPIX2, so the reference stars stay on the same sky (Section 10.3).
- A HISTORY card records the old and new grid size and the padding of each side.

**Frame edges are not data.** The outermost warped pixels of a frame mix its data with the zero fill around it, because the interpolation reaches outside the frame, so their values are attenuated but not 0. The writer does not merge them. On each row and each column of a frame on the grid, it drops a margin of samples inward from the first and the last non-zero pixel: none for nearest-neighbour interpolation, 1 for bilinear, area and Lanczos4, 2 for bicubic, which also overshoots on the next pixel. The margin follows align.interpolation and is not stored. Under CFADRIZ the margin is 0: a drizzled pixel is the weighted mean of the drops that reach it, so a partly covered pixel keeps its full value with a lower weight. An edge at the border of the grid is kept, since the frame may continue beyond it. The footprint that decides the expansion leaves the same samples out, so a grid never grows for pixels whose data is then dropped. Without this rule, an edge value that a grid expansion moves inside the image would be kept for good by the history test (Section 10.6).

The writer refuses an append that would more than double the grid on either axis in one expansion, since that usually means a wrong registration, unless the user forces it. It also refuses one when the accumulator of the expanded grid would not fit in the memory it may use nor in temporary files (Section 10.10). In both cases the user can crop the batch to the grid instead. In a drizzled file, the grid and the padding of each side stay multiples of grid.oversample, the step being rounded up to one, and the footprint reaches half a frame pixel beyond the frame's outer pixel centres at a scale above 1.

Expansion never shrinks the grid. Cropping is a separate operation, out of scope here.

### 10.6 Rejection

Rejection runs in two stages, in this order:

1. **Batch-local rejection.** Siril's existing rejection is run on the batch, unchanged, with reject.method and reject.params. Siril stops rejecting once 4 or fewer samples remain at a pixel, so a batch of 4 frames or fewer gets no batch-local rejection, and only the history test applies.
2. **History test.** This is the only rejection step that LSF adds. Each surviving sample x of a new frame is compared with the accumulated history, where it holds enough samples (see "A thin history is not tested" below). Let σ\_pred be the predicted spread of a new sample about MEAN (Section 10.7). The sample is tested with

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

**A thin history is not tested.** The history test runs only where the history holds at least 5 effective samples, n\_eff = W^2 / W2 (Section 10.7), counted before the batch is merged. This is Siril's own limit: its rejection stops once 4 or fewer samples remain. Elsewhere the samples that survive the batch-local rejection are merged as they are. A few samples cannot tell an outlier from a bad MEAN. When they agree on a wrong value, M2 stays small and every later good sample would be rejected, so the pixel would never heal: edge values that a grid expansion moved inside the image, an outlier in a short first batch, hot pixels. With weights, one heavy sample among light ones counts as about one, however many there are. The limit is a constant of the profile and is not stored. A first batch of fewer than 5 frames therefore leaves its next append untested as well, where nothing else covers its pixels.

### 10.7 Derived statistics

Siril's existing stacking never computes a weighted per-pixel variance. Its rejection uses unweighted statistics of the pixel stack, and the weights enter only the final mean. The statistics below serve only the history test (Section 10.6) and the exported noise and confidence maps.

With W, M2 and W2 the LSF\_ACC planes and N the LSF\_COUNT value at a pixel and channel, the estimator depends on weight.scheme.

**Reliability weights (none, wfwhm, nbstars).** Here the weights express preference, not noise, and the samples are assumed to share a common variance:

- Effective sample size: n\_eff = W^2 / W2.- Unbiased variance of one sample: s^2 = M2 / (W − W2 / W). Undefined with fewer than two contributing samples.
- Variance of the mean: s^2 · W2 / W^2.
- Predicted variance of a new sample about the mean: σ\_pred² = s^2 · (1 + W2 / W^2).

**Precision weights (noise).** Here each weight is proportional to 1 / σ² of its frame, and the weights carry the noise model:

- Unit-weight variance: σ0^2 = M2 / (N − 1). Undefined when N < 2.
- Variance of the mean: σ0^2 / W.
- Predicted variance of a new sample with weight w\_new (the per-sample weight under CFADRIZ, Section 10.3): σ\_pred² = σ0^2 · (1 / w\_new + 1 / W).
- n\_eff = W^2 / W2 is still reported as an effective-integration map.

The two sets agree when all weights are equal. The weight constants are frozen (Section 10.3), so the weights keep the same scale across appends and these estimators stay valid as the stack grows.

### 10.8 Exports

- **Flattened master** (lsf\_export). Read HDU 0, remove the LSF keywords and save it as Siril saves any image. The result is a plain FITS file. If output normalization is enabled, the primary is already the normalized image. The export keeps a HISTORY card naming the source LSFUID and LSFGEN as provenance.
- **Noise and confidence maps** (lsf\_maps). Derived from the accumulator with the statistics of Section 10.7, never stored, and saved as ordinary FITS images with the primary's header (LSF keywords removed) and a HISTORY card naming the source LSFUID and LSFGEN:
  - **neff**: the effective number of frames, W^2 / W2;
  - **noise**: the standard deviation of the mean;
  - **sigma**: the standard deviation of one frame (of a unit-weight frame under noise weighting);
  - **snr**: MEAN divided by the standard deviation of the mean.

  noise and sigma are in the units of the Stack: divided by LSFNMAX - LSFNMIN when LSFNORM = T. A pixel where a value is undefined (no data, or too few samples) is 0. If LSF\_REJ is present, lsf\_maps -rejlo and -rejhi export the low and high rejection maps the same way, as the uint16 counts LSF\_REJ stores, with a HISTORY card giving LSFRGEN. Unlike Siril's own rejection maps, they are not divided by a number of frames: the samples a pixel received since LSFRGEN are not stored.

  lsf\_maps -exp exports the exposure per pixel, in seconds, with one channel per channel of the Stack. It is the exposure map when the file has one. Otherwise, when every ledger row has the same known exposure, it is COUNT times that exposure, and the HISTORY card says so. Otherwise the exposure per pixel is not known, and the export is refused.

The exports are not LSF files, so Siril's FITS preferences apply to them like to any image it saves: its file extension, and its FITS compression when enabled (the file is then named with .fz added, and floating-point data is quantized unless the quantization is 0). With compression off, the flattened master is bit-identical to HDU 0.

### 10.9 Undo

The format keeps no per-batch history planes, because they would grow the file with every append. To be able to undo an append, copy the file before appending. To drop a bad night, keep the copy made before it and continue appending to that copy.

### 10.10 Implementation notes

- **Atomic replace on Windows.** A plain rename fails when the target exists, so it needs MoveFileEx with MOVEFILE\_REPLACE\_EXISTING. It still fails while another program, such as a viewer, holds the file open. A file that was just written can also be held for a moment by a virus scanner, so a sharing or access error is retried for a short time first. If it persists, the writer must report this and keep the temporary file, not lose the commit.
- **Accumulators larger than memory.** An accumulator that would take more than three quarters of the memory the writer may use, counting its optional planes, is kept in memory-mapped temporary files in Siril's swap directory, one per plane, which the system removes when they are closed, even after a crash. The system keeps in memory the parts in use, so the merge, the writer and the expansion work as they do in memory, and a file comes out the same but for rounding in M2. It needs that much free disk during the commit, and is slower once the system runs out of memory to cache the planes. The Stack and the frames being merged stay in the writer's memory. On Linux, a swap directory on a tmpfs defeats this.
- **Long ledger strings and CFITSIO.** CFITSIO 4.7.0 truncates a variable-length string (1PA) that crosses two 2880-byte block boundaries when it is read as text (fits\_read\_col with TSTRING). It returns only the part before the second boundary. ORIGHDR is several kilobytes per row, so some rows are affected in every file, including files written by other tools. A reader that needs the full header reads the column as bytes (TBYTE), using the length from the descriptor (fits\_read\_descript). fitsverify reads these strings the same way and reports "contains non-ASCII text" for the affected rows, because it scans stale bytes in its buffer. That error is a false positive when the stored bytes are printable ASCII. lsf\_compare.py checks the stored bytes directly before it accepts this error. CFITSIO also refuses to write or read as text a string of 28800 characters or more ("ASCII string column is too wide"), which ORIGHDR reaches with a few hundred header cards, as frames registered by Siril often carry: the writer stores variable-length strings as bytes (fits\_write\_col with TBYTE), which gives the same bytes in the file. The transforms of an astrometric registration are computed by the writer and kept by Siril's apply (registration\_args.keep\_transforms), so a stale external reference left in a .seq file does not change what Siril's own seqapplyreg does.

### 10.11 Mosaics

A mosaic is a grid that grows panel by panel. The format needs nothing more than grid expansion (Section 10.5) and the exposure map (Section 3), which gives the integration time of every pixel. Three append options and one choice of weighting make it work.

- **Workflow, panel by panel.** Plate-solve the first panel, then initialize from it, or import its WCS later (lsf\_util -wcs). Append every other panel with lsf\_append -wcsreg when it shares too few stars with the Stack (Section 10.4), with -overlap\_norm and -feather=, and with -expand\_force when a panel more than doubles the grid on one axis. lsf\_maps -exp and -neff show the depth of every pixel.
- **Workflow, all panels at once.** Plate-solve every frame of every panel with registration (seqplatesolve without -noreg), then initialize from that sequence with -overlap\_norm and -feather= (Section 10.1). Later frames are plate-solved and appended with -wcsreg.
- **Weighting.** Use noise weighting. The wFWHM and star counts are measured against LSF\_WREF, the stars of the first batch's reference frame (Section 10.3). A panel away from it covers few or none of them, so its wFWHM is unreliable or undefined, and its star count follows the field rather than the quality. The writer warns about frames that cover fewer than 10 reference stars. A frame under wfwhm that covers none is not merged, with a message saying so. Under nbstars, a batch that covers less than half of the reference frame on average gets a warning.
- **Normalization on the overlap** (lsf\_append -overlap\_norm). The frozen baseline compares a frame's global statistics with those of the first batch's reference frame, which holds for frames of one field. A panel of another field has other global statistics, and its level is matched wrongly, leaving a seam.
  - With -overlap\_norm, the frame of the batch that overlaps the Stack most is matched to it on their overlap: its data pixels where W > 0, against MEAN on the same pixels.
  - Only the level comes from the overlap: the offset, the median of MEAN minus the scaled frame, or for multiplicative normalization the ratio of their medians. The scale stays the frozen one. The scale of Siril's normalization is a dispersion, mostly noise. It compares two frames, but not a frame with the Stack, whose noise is averaged away. A fit of the frame against MEAN would be biased by their different star profiles.
  - The other frames of the batch are normalized onto that frame by their global statistics, as usual, so a batch should show one panel. Each channel is matched to its own, the Stack being equalized already.
  - An overlap under 10000 pixels or 1 percent of the frame falls back to the frozen baseline with a warning. NORM\_MTH records '+ovl' (Section 9).
  - In a test (two panels with 20 percent overlap, the second at 1.15 × level + 0.01), the level of the second panel was off by 9.3σ of the sky without the option, and by 1.8σ with it.
- **Normalization on the overlaps at initialization** (lsf\_init -overlap\_norm), for a first batch of several panels.
  - The frames are first normalized by their global statistics against the reference frame, as usual. Their levels are then matched where they overlap, all frames together, with Siril's overlap solver: each pair of frames gives the median of each on the bins they share, and the solver finds the level of every frame connected to the reference frame by overlaps, which keeps its own.
  - The matching works on bins of 16 by 16 pixels of the grid (larger when the binned frames would take more than a quarter of the memory Siril may use), each the mean of a frame where the bin is all data, in every channel. The result is an offset, or for multiplicative normalization a factor. Only the level is matched; the scale stays the frame's own, for the reasons above. Frames that are not merged, because their weight is 0, take no part.
  - A pair of frames sharing fewer bins than 10000 pixels or 1 percent of either frame's is not used. A frame the remaining pairs do not connect to the reference frame keeps its global normalization, with a warning. The frozen baseline stays that of the reference frame. NORM\_MTH records '+ovl'.
- **Feathering** (lsf\_init and lsf\_append -feather=N, 0 to 2000).
  - It uses Siril's blending masks: a distance to each frame's edge, computed at a tenth of the size, ramped over N (in the units Siril's stack uses: pixels of the masks, so N is about 10 N pixels of the frames).
  - The ramp multiplies the sample's weight (Section 10.3). Siril upscales each block of rows separately, so its weights depend on how a stack is cut in blocks. LSF interpolates each row from the whole mask, which is what Siril gives when a channel is one block.
  - Only the incoming frames can be ramped. The edges of data already in the Stack cannot be feathered after the fact, because only a common factor can be applied to the past (Section 10.3). Where a new panel reaches beyond the old coverage, the old side keeps its hard edge. With -overlap\_norm, what remains there is a step in noise, not in level.
  - FEATHER records the distance (Section 9). The masks of a batch the command registered are removed with it.

## 11. Conformance

- **Any FITS reader.** No requirements. It benefits from C1 to C4.
- **LSF reader.** Locates HDUs by name and cube planes by LSFPLn, validates per Section 8, reads the ledger, supports GZIP\_1, GZIP\_2 and RICE\_1 tile compression (refusing a lossy LSF\_ACC, Section 5.3), and tolerates missing optional keywords, columns and HDUs (LSF\_REJ, LSF\_STARS).
- **LSF writer.** Everything an LSF reader does, plus the pre-merge validation, the commit protocol, preservation of unknown content, and honest maintenance of the primary's standard keywords (Section 7). In Siril, writing is done only by the LSF commands (Section 1).

## 12. Possible Future Improvements

Possible future improvements are kept outside this profile, in `Future work.md` (beside this file), until they are taken up. They include reference stars for mosaic panels (LSF\_WREF beyond the first batch), export to FLIS, nbstack weighting, in-place append and frame filters in the commands.

## 13. Known Limitations

- **Append is Siril-only.** Other tools can open the file and use the primary image, but they are not expected to preserve the extension HDUs, the LSF keywords or the ledger when they save.
- **A foreign save can end the stack's life.** A tool that re-saves the file may drop the extensions or strip the LSF keywords. The result is still a valid image, but the accumulator can no longer be merged into, and continuing the stack means starting again from the original frames. A tool that keeps the extensions but changes the primary is detected and repairable (Section 8).
- **Process exports, not the LSF file.** Edits made to the primary are either undone by repair, or adopted without adjusting M2.
- **Per-pixel integration time needs known exposures.** It is exact with the exposure map, or COUNT times the exposure while every frame has the same one (Section 10.8). A file whose ledger holds a frame of unknown exposure, or that mixed exposures before exposure maps existed, cannot get one. A frame of unknown exposure merged into a file with a map adds 0 s to it (Section 10.2).
- **No built-in undo.** Undo is done by keeping copies (Section 10.9).
- **Drizzle limits.** The scale is a whole number from 1 to 3. Without output normalization, a sky above 1/s² clips the whole Stack at scale s: no star is found in it, so frames can no longer be registered against it. MEAN is not affected, but a file cannot change to output normalization later. A master flat is not used as drizzle weights: they are geometry only. Feathering, -overlap\_norm, -noreg, -wcsreg, distortion correction and a registration through plate solving are refused on a drizzled file, as they take the holes of its frames for dark pixels or are not tested with them.
