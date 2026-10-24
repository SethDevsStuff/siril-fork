# LSF UI integration plan (draft, 2026-10-22, reviewed against 8341c3d8e on 2026-10-23)

Scope: GTK4 GUI (src/gui-gtk4) for the LSF commands that exist today: lsf_init, lsf_append, lsf_export, lsf_util (check, repair, adopt, compress, decompress, copy) and lsf_maps. Nothing in src/gui-gtk4 references LSF yet.

Section A covers what the GUI needs from the core first, section B the user-facing work, section C the decisions to make before starting. The code and spec as of 8341c3d8e are the reference for this plan.

## Progress (2026-10-23, uncommitted)

- **B5 done.**
  - `fits_file_is_lsf()` (io/image_format_fits.c); `savefits()` and the stacking output refuse to overwrite an LSF file.
  - Opening an LSF sets `single.is_lsf` and strips the LSF keywords from the view, so Save As writes a plain image. Save is disabled, the suggested name drops `.lsf`, and the file label shows "(Living Stack, read-only)".
  - Tested headless: save and stacking onto the LSF name are refused, and the file stays byte-identical. Not yet checked in the GUI.
- **A1 done.** `lsf_read_info()` / `lsf_info` (lsf_fits.c) share `open_and_validate()` with the full validation, in a light mode: checks 1–4, 6 and 7, no planes or checksums, failure text instead of a log line. It adds the ledger summary and `no_append`. Also exposed as `lsf_util -info`.
- **A2 + A3 done.**
  - `struct lsf_cmd_args` is in lsf.h, with `lsf_cmd_args_new()`, `lsf_start()`, `lsf_cmd_args_to_command()`, `end_func` (run on the main thread after the processing thread is released), `export_out` (the "+ image" export) and `load_result` (the new Stack becomes the result image, the sequence staying loaded).
  - The command parsers use the same constructor; behaviour unchanged (regression-tested).
- **B1 done, waiting for a GUI test.** Result drop-down (Image / LSF / LSF + image) and the `stack_lsf_grid` panel (framing, compression, history test, rejection maps). The rejection maps option is shown whatever the rejection method, since the history test rejects on appends, and its tooltip explains the stored counts.
  - **File name:** a pattern, like the stacking result. The default is `$seqname$$FILTER:%s$`.
    - It's parsed with `update_header_and_parse` (no-fail, creating directories) against the reference frame's header, since the Stack doesn't exist yet.
    - `.lsf.fits` is added: a `.fit`/`.fits`/`.fts` ending and trailing `_`/`-` are removed first.
    - Image mode and LSF mode keep separate patterns. `stack_set_image_output_name()` updates the image one, from a sequence load or the stack-name preference, without touching the LSF one.
    - The command line logs the parsed name, since `lsf_init` takes a literal path. The other modes are hidden as planned. The equivalent command is logged.
  - **A4, partly:** LSF mode uses the sequence's own selection, saved before the run, and hides the stacking filters, as the multipoint stack does. Passing the filtered list is still open. `lsf_process_batch` writes the batch's `.seq` (statistics cache), so a selection changed in memory would leak into the user's file.
- **A5 done** for init, as a confirmation dialog.
- **A6 done for append.** `lsf_append_preflight()` (lsf_commands.c):
  - memory: refused when the accumulator alone exceeds the limit, warned above 75%;
  - disk: the new file beside the old one, plus the registered `lsf_r_` frames in the working directory (one float32 grid-size frame per selected image, unless noreg), summed when both are on one file system.

  The grid growth is still unknown before the registration. The init pre-flight is not done.
- **B2 done, waiting for a GUI test.** Result entry "Append to a Living Stack", with `stack_lsf_append_box`: file entry and browse button, info label, grid grow/crop, step, force, noreg, compression keep/none/gzip2/gzip1, backup (on).
  - **On choosing a file:** `lsf_read_info`, then `fill_widgets_from_lsf()`. Rejection, normalization, weighting and output normalization are locked. The user's stacking preferences are saved and restored around the fill.
  - **History test:** editable, starting from the file's values. A change is confirmed with the decided notice and passed only for the values that changed.
  - **Rejection maps:** the LSF option becomes "Add rejection maps…", shown only when the file has none.
  - **Checks:** `lsf_check_batch()` refuses a channel mismatch and warns when the reference frame's Section 7 keywords differ, through `lsf_frame_mismatch()`, the same check lsf_append makes. Problems show in the info label and disable Run when fatal.
  - **After the run:** the Stack becomes the result image and the file is read again. The last appended file is proposed again for the session only; the per-directory memory is not done.

## A. Core prerequisites

These come first, because every UI item below depends on them.

1. **A light reader for pre-filling the UI.** `lsf_open_and_validate` reads HDU 0, every LSF_ACC/LSF_COUNT/LSF_REJ plane and the checksums. That can be several GB and many seconds, which is far too slow to run every time the user picks a file in a chooser. Add `lsf_read_info(path, lsf_info *out)`, which reads only:
   - the primary keywords;
   - LSF_SETTINGS;
   - the dimensions and ZCMPTYPE of each HDU;
   - a ledger summary: rows, batches, generation, first and last DATE-OBS, total exposure, filters.

   No planes and no checksums. It reports the cheap checks (1, 2, 3, 7) and leaves the expensive ones (5, 8) to the actual run. The summary lsf_util -check prints (`print_summary`) can then be built from this struct as well.

   It must also report the reasons a file cannot be appended to, so that the UI can say so before the run:
   - grid.oversample other than 1 (lsf_append now refuses a drizzled grid);
   - cfa.mode = CFADRIZ (refused by the settings reader);
   - history settings that are not valid (`lsf_settings_from_rows` now refuses thresholds that are not above 0, and a margin below 1).

   lsf_export now runs the full validation as well (it tolerates an edited Stack), so Export reads the whole accumulator and needs a progress bar, like append.
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
6. **Disk and memory pre-flight.**
   - **Disk.** An append writes the new file next to the old one before replacing it. The .bak is the old file itself, renamed in the same step as the replace (`ReplaceFileW` on Windows, a hard link elsewhere), not a copy. So the free space needed is one new file: `acc_bytes + count_bytes + rej_bytes` plus HDU 0 and the grid growth when known. With a backup, that space stays used afterwards. Check it with `test_available_space()` (core/OS_utils.h) and warn when it is short.
   - **Memory.** `lsf_process_batch` now counts the accumulator against Siril's memory limit. A batch fails with "not enough memory" when the accumulator (rx · ry · nchan · (16 B + 2 B, or 6 B with LSF_REJ)) leaves too little for the frame rows. The GUI can compute this from `lsf_read_info` before the run, including the expanded grid when it is known. It should say what to change (memory limit, or framing = current) instead of letting the job fail.
   - **Failed replace.** `lsf_write_file` now reports two cases. Both should be dialogs, not only log lines:
     - The file could not be replaced: it is unchanged, and the new file is kept under the temporary name (Section 10.10).
     - The backup rename succeeded but the replace did not, and the rename back also failed: the file now exists only as the .bak, and the user must rename it back.

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
  | `stack_lsf_rejmaps` (its own option, not the stacking `rejmaps_checkbutton`, which is hidden with the merge option) | -rejmaps |
  | `entryresultfile` | file path, forced to end in .lsf.fits |
  | `checkbutoverwrite` | -overwrite |

- **Framing.** Replace `check_maximize_framing` with a drop-down: max (default), current, min, cog. Unlike stacking, lsf_init applies the registration itself, so:
  - framing is meaningful even when the registration is not yet applied;
  - the "Registration data found, apply first" confirmation in `start_stacking()` must be skipped in LSF mode.
- **Hidden in LSF mode:** feather/blend, overlap normalization, upscale at stacking (until CFADRIZ), force 32-bit (the LSF is always float32), merge low/high rejmaps (LSF_REJ keeps low and high), force normalization, and the MPP pages.
- **New LSF-only panel** (a small expander, visible only in LSF modes):
  - history test: low, high and margin, with an "Same as batch sigma" toggle that is on by default. As lsf_init now does, the toggle falls back to 3.0 when the batch method has no sigma threshold or its sigma is 0. The spin buttons need limits matching the parser: thresholds above 0 (0 would reject every sample on that side), margin at least 1;
  - compression: none / gzip2 / gzip1.
- **"+ image".** After the commit, run `lsf_export_master` to `<name>.fit` (Siril's extension and compression preferences apply, Section 10.8). In both LSF modes, load the Stack into the display afterwards, as stacking does with its result. The display should be marked as an LSF view (B5).
- **Default name.** `$seqname$_<FILTER>.lsf.fits`. The file is one per filter, so putting the filter in the name helps the user pick the right file later.

### B2. Appending to an LSF

The registration does not need a separate step, because the append already does it. `register_against_stack()` runs Siril's global star alignment with `external_ref_path = <lsf>` (HDU 0), using the findstar, transform, interpolation and min-pairs settings stored in the file. It then uses the footprint to plan the grid expansion and applies the registration onto the expanded canvas. The apply step is tied to the expansion, so it cannot move to the Registration tab.

Recommended flow: a fourth result entry on the stacking tab, **Append to LSF…**, with a file chooser. On choosing a file:

1. Run `lsf_read_info` (A1) and show an info box: frames, batches, generation, grid size and scale (grid.oversample), channels, filter, total exposure, last session date, compression, rejection maps. When the file cannot be appended to by this version (a drizzled grid, CFADRIZ, invalid history settings), say why in the box and disable Run.
2. **Pre-fill and lock the widgets.**
   - **Facts** are filled and made insensitive: normalization method, estimator, RGB equal, weighting, output normalization. They are frozen in the file, and changing them would make the accumulator meaningless.

     Output normalization is read from the file each time, not assumed from init. A commit that finds no valid dynamic range now turns it off for good: `lsf_build_primary` sets `outnorm` to FALSE, and later appends do not turn it back on.
   - **Defaults** are filled from the file: the rejection method and sigmas, and the history low, high and margin.
     - The history values can be overridden on append (`-hist_*`). An override is saved in the file and becomes the default of later appends (decided 2026-10-23, implemented in 8341c3d8e).
     - When the user changes one of them, show a notice before the run: the change cannot affect the frames already in the stack, and it is saved for the appends that follow. lsf_append logs the same thing.
     - The rejection method and sigmas cannot be overridden yet (C2). If they become overridable, the same rule applies.
     - The change is only logged, so the file itself does not keep which thresholds earlier batches used. lsf_append should also add a HISTORY card to HDU 0 ("LSF gen N: history test low/high/margin a/b/c, was x/y/z"), as it already does for grid expansion, so the LSF Manager (B4) can show it.
3. **Compatibility pre-check** against the loaded sequence: channel count, sensor class and Bayer pattern, binning, filter, and pixel size and focal length (Section 7). Show the problems inline and disable Run on a hard mismatch. Today the user would only see them in the log after the job starts.
4. **Append-only widgets:**
   - framing: max (grow the grid) or current (crop to it), since min and cog do not apply on append;
   - expansion step, default 64;
   - "allow growing more than 2×" (-expand_force);
   - "frames already on the grid" (-noreg);
   - "add rejection maps" (-rejmaps), shown only when the file has none, with a note that the counts start from this append;
   - compression: keep (default), change, or remove;
   - **keep a backup (.bak), on by default in the GUI**, because there is no undo (Section 10.9).
   - Later: registration by stars or by WCS. Section 10.4 of the spec describes `lsf_append -wcsreg`, which is not implemented yet.
5. Run, then reload the Stack in the display and refresh the info box.

The tab should also remember the last LSF used for a sequence directory, so a second night is a single click.

### B3. Registration tab: register against an LSF (optional, phase 2)

The Registration tab already supports an external reference for global and 2-pass star alignment (`reg_reference_checkbutton`, `reg_referencefilechooser_box`, `regargs->external_ref_path`). Choosing an .lsf.fits there already works mechanically, because HDU 0 loads as an image. What the tab adds is inspection: the user sees the registration results in the plot and the frame list, and can drop bad frames before an append that cannot be undone.

To make this a real flow:

- When the external reference is an LSF file, pre-fill transform, interpolation and min pairs from align.*.
- Use the star detection from align.findstar for this run only, through `registration_args.sf_params` (added in 8341c3d8e). The append already does this; the user's preferences are not touched.
- Force "compute only" (no output sequence). The apply step must stay in the append, for the expansion canvas.
- Store, with the regdata, a tag naming what it was computed against: LSFUID, LSFGEN and the grid size.
- The append itself now leaves the batch's .seq file and its cached star lists as the user made them (`no_starlist`, no `writeseqfile`). This flow would be the only way regdata against an LSF ends up in the user's .seq, which is one more reason for the tag:
  - a later ordinary registration replaces that regdata;
  - stacking it as a normal sequence would be wrong.

  Use `no_starlist` here too, so the star lists detected with the file's settings do not replace the user's cached ones.
- lsf_append gains a `-usereg` mode, and the GUI's Append offers "Use existing registration against this file", which skips pass 1. It is refused when the tag does not match the file (another file, or a later generation whose grid has grown).

Framing does not belong in this tab, because it is decided when the registration is applied, which the append does (B2).

Recommendation: ship B2 without this, then add B3 if inspecting before committing proves useful.

### B4. Tools > Living Stack

The user can pick each action from a menu, but every action needs the same file chooser and the same summary. Recommended: one **LSF Manager** dialog, opened from a Tools entry in `tools_menu` and from the stacking tab's info box. Its parts:

- **Header:** a file chooser and the `lsf_read_info` summary.
- **Ledger view:** a table with one row per frame (batch, date, filter, exposure, FWHM, weight, noise, stars) and a per-batch grouping, so the user can see which night contributed what. Sorting on FWHM or weight helps to spot a bad session.
  - Under wfwhm and nbstars weighting, a frame at or beyond the first batch's worst is no longer dropped. It gets the weight floor (`LSF_WEIGHT_FLOOR`, 0.01 before normalization, Section 10.3). Mark rows at the floor: a batch with many of them was worse than anything in the first batch, which the user may not have noticed.
  - Show the HISTORY cards of HDU 0 (expansions, repairs, adoptions, history-test changes) as a per-generation log beside the table.
- **Check:** the full validation (lsf_util -check), with the result shown per check, not only in the log.
- **Repair / Adopt:** enabled only when check 8 failed. Each needs a one-line explanation:
  - Repair discards the outside edit.
  - Adopt keeps the edit as the new mean, without adjusting M2. NaN or infinite edited pixels are not adopted and keep the accumulator's value; adopt now counts them. Show the count in the result, not only the log.
- **Compress / Decompress:** compression drop-down (gzip2, gzip1), showing the current codec and stored size of each HDU, which are already in `lsf_file`.
- **Save copy as…:** -out, which with no other action makes a copy. Useful as a manual undo point before a risky append.
- **Export master…:** lsf_export. It now exports a Stack edited outside Siril with the edit, and warns. When check 8 fails, the dialog should say so next to the button, and offer Repair or Adopt first.
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
- **Cancellation:** the commit is atomic, so stopping before the rename is safe. Check that the merge and write loops test `get_thread_run()`, and that a stop removes the temporary file and the `lsf_r_` sequence. The batch's own .seq is no longer written during an append, so those are the only files a stop can leave behind.
- **Plate-solving the Stack (needed by -wcsreg, Section 10.4).** The astrometry dialog on an LSF view should offer "write WCS into the LSF" instead of a save, through a new `lsf_util -setwcs`-style path. Otherwise WCS registration on append can only ever use a WCS inherited from the first batch.
- **CFADRIZ (Sections 10.3 and 10.4).** Today the reader refuses cfa.mode = CFADRIZ, and lsf_append refuses grid.oversample other than 1. When CFADRIZ lands, the LSF panel needs the drizzle scale, pixfrac and kernel. Leave room for them in the panel layout rather than re-using the hidden "upscale at stacking" widgets.
- **Live stacking (later).** livestacking_gui is the obvious long-term consumer: commit to an LSF every N frames. Keep the A2 API free of tab-specific state so it can be driven from there.

## C. Decisions needed

1. **Where the LSF modes live on the stacking tab.**
   - (a) A result drop-down (Image / LSF / LSF + image / Append), as in B1 and B2. Recommended.
   - (b) A separate LSF expander that switches the tab into LSF mode.
2. **Rejection override on append.** The rule for overrides is decided: an override is saved as the new default (B2). The history thresholds follow it. What is still open is whether the rejection method and sigmas can be overridden too; lsf_append has no option for them.
   - (a) Lock the widgets, matching the code today.
   - (b) Add `rej`/sigma override to lsf_append, so the widgets stay editable. The ledger, or a HISTORY card as proposed for the history test in B2, would then need to record per batch the rejection that was used.
3. **Registration tab flow (B3):** phase 2, or not at all.
4. **Tools:** one LSF Manager dialog (recommended), or separate menu entries per action.
5. **Maps:** in-memory viewing (needs the refactor in B4), or saving to files only.

Suggested order: B5 → A1 → A2 + A3 → B1 → A4 + A5 + A6 → B2 → B4 → B6 → B3.
