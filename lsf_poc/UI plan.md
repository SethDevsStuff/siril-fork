# LSF UI integration plan (draft, 2026-10-22)

Scope: GTK4 GUI (src/gui-gtk4) for the LSF commands that exist today: lsf_init, lsf_append, lsf_export, lsf_util (check, repair, adopt, compress, decompress, copy) and lsf_maps. Nothing in src/gui-gtk4 references LSF yet.

Section A covers what the GUI needs from the core first, section B the user-facing work, section C the decisions to make before starting.

## A. Core prerequisites

These come first, because every UI item below depends on them.

1. **A light reader for pre-filling the UI.** `lsf_open_and_validate` reads HDU 0, every LSF_ACC/LSF_COUNT/LSF_REJ plane and the checksums. That can be several GB and many seconds, which is far too slow to run every time the user picks a file in a chooser. Add `lsf_read_info(path, lsf_info *out)`, which reads only:
   - the primary keywords;
   - LSF_SETTINGS;
   - the dimensions and ZCMPTYPE of each HDU;
   - a ledger summary: rows, batches, generation, first and last DATE-OBS, total exposure, filters.

   No planes and no checksums. It reports the cheap checks (1, 2, 3, 7) and leaves the expensive ones (5, 8) to the actual run. The summary lsf_util -check prints (`print_summary`) can then be built from this struct as well.
2. **A public entry point the GUI can call.** `struct lsf_cmd_args` and the `*_run` functions are private to lsf_commands.c, and the only way in is the word[] parser. Two options:
   - Expose `lsf_init_params` / `lsf_append_params` / `lsf_util_params` plus `lsf_start_*()` in lsf.h, and have the `process_lsf_*` parsers fill them. The GUI fills the same structs.
   - Or have the GUI build a command line and run it through processcommand. This is cheaper, but it loses typed errors and the completion hook.

   Recommended: the first. The GUI should still log the equivalent command line, so a GUI run can be repeated in a script.
3. **A completion hook.** `lsf_worker` ends with `siril_add_idle(end_generic, NULL)`. The GUI needs its own idle callback with the result, to load the new Stack into the display, refresh the LSF info panel and report "frames merged / dropped / grid grown". Add an optional `idle_func` + user data to the params.
4. **Frame selection.** `load_batch()` loads the sequence from disk, so lsf_init and lsf_append use the selection saved in the .seq file, not the one in memory.
   - The stacking tab's image filters (FWHM, roundness, quality, % / k-sigma) build their own include list in `stackparam`. LSF currently ignores them.
   - Add an optional `int *image_indices, nb` to the batch params, used instead of `seq->imgparam[i].incl`. The GUI fills it from the same filter code as stacking.
   - At the very least, `writeseqfile(&com.seq)` before starting, so manual selection changes are honoured.
5. **Low-stack warning as a question.** Section 10.1 asks the writer to warn when the first batch has fewer than 25 frames. Today this is a log line. In the GUI it should be a `siril_confirm_dialog` before the run, so the check belongs in a pre-flight function that the GUI calls, not in the worker.
6. **Disk-space pre-flight.** An append writes a full temporary copy, plus the .bak when asked. Estimate the size from `acc_bytes + count_bytes + rej_bytes` (and the grid growth when known) and warn when the volume's free space is short. On Windows, a rename that keeps failing leaves the temporary file (Section 10.10). The GUI should show that as a dialog, not only as a log line.

## B. UI work

### B1. Stacking tab: create an LSF

Add a result selector to the output expander (`stack_expander_output`): **Image** (today's behaviour), **LSF file**, **LSF file + image**. A drop-down is the better fit, because B2 adds a fourth entry.

When an LSF entry is selected:

- **Method.** Locked to "Average stacking with rejection" (`comboboxstack_methods` insensitive). lsf_init supports only that method.
- **Widgets that map to lsf_init:**

  | Widget | lsf_init |
  |---|---|
  | `comborejection` + sigma spins | rej method, sig low/high |
  | `combonormalize` | -norm= / -nonorm |
  | `checkfastnorm` | -fastnorm |
  | `check_RGBequal` | -rgb_equal |
  | `comboweighing` | -weight= (hide nbstack, which is not supported) |
  | `check_normalise_to_max` | -output_norm |
  | `rejmaps_checkbutton` | -rejmaps |
  | `entryresultfile` | file path, forced to end in .lsf.fits |
  | `checkbutoverwrite` | -overwrite |

- **Framing.** Replace `check_maximize_framing` with a drop-down: max (default), current, min, cog. Unlike stacking, lsf_init applies the registration itself, so:
  - framing is meaningful even when the registration is not yet applied;
  - the "Registration data found, apply first" confirmation in `start_stacking()` must be skipped in LSF mode.
- **Hidden in LSF mode:** feather/blend, overlap normalization, upscale at stacking (until CFADRIZ), force 32-bit (the LSF is always float32), merge low/high rejmaps (LSF_REJ keeps low and high), force normalization, and the MPP pages.
- **New LSF-only panel** (a small expander, visible only in LSF modes):
  - history test: low, high and margin, with an "Same as batch sigma" toggle that is on by default;
  - compression: none / gzip2 / gzip1.
- **"+ image".** After the commit, run `lsf_export_master` to `<name>.fit` (Siril's extension and compression preferences apply, Section 10.8). In both LSF modes, load the Stack into the display afterwards, as stacking does with its result. The display should be marked as an LSF view (B5).
- **Default name.** `$seqname$_<FILTER>.lsf.fits`. The file is one per filter, so putting the filter in the name helps the user pick the right file later.

### B2. Appending to an LSF

The registration does not need a separate step, because the append already does it. `register_against_stack()` runs Siril's global star alignment with `external_ref_path = <lsf>` (HDU 0), using the findstar, transform, interpolation and min-pairs settings stored in the file. It then uses the footprint to plan the grid expansion and applies the registration onto the expanded canvas. The apply step is tied to the expansion, so it cannot move to the Registration tab.

Recommended flow: a fourth result entry on the stacking tab, **Append to LSF…**, with a file chooser. On choosing a file:

1. Run `lsf_read_info` (A1) and show an info box: frames, batches, generation, grid size, channels, filter, total exposure, last session date, compression, rejection maps.
2. **Pre-fill and lock the widgets.**
   - **Facts** are filled and made insensitive: normalization method, estimator, RGB equal, weighting, output normalization. They are frozen in the file, and changing them would make the accumulator meaningless.
   - **Defaults** are filled from the file: the rejection method and sigmas, and the history low, high and margin. The history values can already be overridden on append (`-hist_*`). Rejection cannot (see C2). An override is saved in the file and becomes the default of later appends (decided 2026-10-23). When the user changes one of these values, show a notice before the run: the change cannot affect the frames already in the stack, and it is saved for the appends that follow. The command logs the same thing.
3. **Compatibility pre-check** against the loaded sequence: channel count, sensor class and Bayer pattern, binning, filter, and pixel size and focal length (Section 7). Show the problems inline and disable Run on a hard mismatch. Today the user would only see them in the log after the job starts.
4. **Append-only widgets:**
   - framing: max (grow the grid) or current (crop to it), since min and cog do not apply on append;
   - expansion step, default 64;
   - "allow growing more than 2×" (-expand_force);
   - "frames already on the grid" (-noreg);
   - "add rejection maps" (-rejmaps), shown only when the file has none, with a note that the counts start from this append;
   - compression: keep (default), change, or remove;
   - **keep a backup (.bak), on by default in the GUI**, because there is no undo (Section 10.9).
   - Later, from M5: registration by stars or by WCS.
5. Run, then reload the Stack in the display and refresh the info box.

The tab should also remember the last LSF used for a sequence directory, so a second night is a single click.

### B3. Registration tab: register against an LSF (optional, phase 2)

The Registration tab already supports an external reference for global and 2-pass star alignment (`reg_reference_checkbutton`, `reg_referencefilechooser_box`, `regargs->external_ref_path`). Choosing an .lsf.fits there already works mechanically, because HDU 0 loads as an image. What the tab adds is inspection: the user sees the registration results in the plot and the frame list, and can drop bad frames before an append that cannot be undone.

To make this a real flow:

- When the external reference is an LSF file, pre-fill transform, interpolation and min pairs from align.*. Pre-fill the star detection from align.findstar for this run only, without writing it to the user's preferences (`register_against_stack` already saves and restores `com.pref.starfinder_conf`).
- Force "compute only" (no output sequence). The apply step must stay in the append, for the expansion canvas.
- Store, with the regdata, a tag naming what it was computed against: LSFUID, LSFGEN and the grid size.
- lsf_append gains a `-usereg` mode, and the GUI's Append offers "Use existing registration against this file", which skips pass 1. It is refused when the tag does not match the file (another file, or a later generation whose grid has grown).

Framing does not belong in this tab, because it is decided when the registration is applied, which the append does (B2).

Recommendation: ship B2 without this, then add B3 if inspecting before committing proves useful.

### B4. Tools > Living Stack

The user can pick each action from a menu, but every action needs the same file chooser and the same summary. Recommended: one **LSF Manager** dialog, opened from a Tools entry in `tools_menu` and from the stacking tab's info box. Its parts:

- **Header:** a file chooser and the `lsf_read_info` summary.
- **Ledger view:** a table with one row per frame (batch, date, filter, exposure, FWHM, weight, noise, stars) and a per-batch grouping, so the user can see which night contributed what. Sorting on FWHM or weight helps to spot a bad session.
- **Check:** the full validation (lsf_util -check), with the result shown per check, not only in the log.
- **Repair / Adopt:** enabled only when check 8 failed. Each needs a one-line explanation:
  - Repair discards the outside edit.
  - Adopt keeps the edit as the new mean, without adjusting M2.
- **Compress / Decompress:** compression drop-down (gzip2, gzip1), showing the current codec and stored size of each HDU, which are already in `lsf_file`.
- **Save copy as…:** -out, which with no other action makes a copy. Useful as a manual undo point before a risky append.
- **Export master…:** lsf_export.
- **Maps:** check boxes for neff, noise, sigma, snr, rejlo and rejhi (the last two only with LSF_REJ), plus an output prefix to save them.
- **View plane:** loads one map into the display as a new unsaved image, without writing files. It needs a small refactor of `lsf_maps_run` into a compute step (into a `fits*`) and a save step. To keep memory low, it should read only the planes that map needs (MEAN, W, M2, W2, COUNT), not run the full validation.

  Raw W and COUNT views are also worth offering: they show where the coverage is thin at the edges after a grid expansion.

### B5. Protecting LSF files opened as images (do first)

File > Open of an .lsf.fits loads HDU 0 as a normal image. A plain Save then calls `savefits(com.uniq->filename, gfit)` (save_dialog.c), which rewrites the file without its extensions. Siril would then be the "foreign save" that ends the stack's life (Section 13), and it takes a single Ctrl+S.

- On open, detect `LSF = T` and set a flag on the loaded image.
- While the flag is set:
  - redirect Save to Save As, suggesting `<name>.fit`, and refuse a target that is an LSF file;
  - show an info bar: "Living Stack: read-only view. Use Tools > Living Stack to export or maintain it."
- Add a "Living Stack (*.lsf.fits)" filter to filefilters.ui, for the choosers in B1 to B4.

This needs no other LSF UI and protects users today, since they can already open these files.

### B6. Smaller items

- **Preferences (LSF page or a section of Stacking):** default compression, backup on append, expansion step, history margin, the low-stack warning threshold (25), and the default framing.
- **Cancellation:** the commit is atomic, so stopping before the rename is safe. Check that the merge and write loops test `get_thread_run()`, and that a stop removes the temporary file and the `lsf_r_` sequence.
- **Plate-solving the Stack (ties to M5 B5).** The astrometry dialog on an LSF view should offer "write WCS into the LSF" instead of a save, through a new `lsf_util -setwcs`-style path. Otherwise WCS registration on append can only ever use a WCS inherited from the first batch.
- **CFADRIZ (M5 B4).** When it lands, the LSF panel needs the drizzle scale, pixfrac and kernel. Leave room for them in the panel layout rather than re-using the hidden "upscale at stacking" widgets.
- **Live stacking (later).** livestacking_gui is the obvious long-term consumer: commit to an LSF every N frames. Keep the A2 API free of tab-specific state so it can be driven from there.

## C. Decisions needed

1. **Where the LSF modes live on the stacking tab.**
   - (a) A result drop-down (Image / LSF / LSF + image / Append), as in B1 and B2. Recommended.
   - (b) A separate LSF expander that switches the tab into LSF mode.
2. **Rejection override on append.** Section 6 says reject.* defaults are recalled "unless the user overrides them for that append", but lsf_append has no option for it.
   - (a) Lock the widgets, matching the code today.
   - (b) Add `rej`/sigma override to lsf_append, so the widgets stay editable. The ledger would then need a per-batch record of the rejection that was used.
3. **Registration tab flow (B3):** phase 2, or not at all.
4. **Tools:** one LSF Manager dialog (recommended), or separate menu entries per action.
5. **Maps:** in-memory viewing (needs the refactor in B4), or saving to files only.

Suggested order: B5 → A1 → A2 + A3 → B1 → A4 + A5 + A6 → B2 → B4 → B6 → B3.
