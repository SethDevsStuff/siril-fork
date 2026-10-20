#!/usr/bin/env python3
"""
lsf_export.py -- export a Living Stack Format (LSF) file's flattened master
(and, optionally, a per-pixel noise/variance map) into a plain FITS
file that ordinary tools -- ones with no idea LSF exists -- can open.

Author: (fill in your name here before distributing)
License: (fill in a license here before distributing, e.g. GPL-3.0)
Contact: (fill in a contact / forum / repo URL here before distributing)

This corresponds to the "Export a flattened master" and "Export a
confidence/noise map" operations in live_stack.md, Section 9:

    Export a flattened master: mean = sum / weight. This is an ordinary
    float32 image, identical in precision to what any conventional
    stacking tool would export -- not a lossy approximation, just a
    division.

    Export a confidence/noise map: derived directly from m2, weight and
    the sum of squared weights (LSF_W2) -- a per-pixel variance or
    effective-integration map.

Since LSF already stores the flattened mean directly as the "Stack" image
plane (see lsf_stack.py), exporting the master is actually just a format
conversion of that one plane; the noise map is where sumsq/m2 and weight
get used.

The input is an LSF FITS file (live_stack_fits.md). The export
is a PLAIN FITS image: every LSF identity keyword (LSF, LSFVER, LSFUID, LSFGEN, LSFNORM, ...) is
removed so no tool mistakes it for an LSF file, and a HISTORY card records which LSF file and
generation it came from. Exporting never validates the input: the Stack must stay readable
even when the file's integrity binding is broken (C1).

WHICH MAP? Two different noise things are derivable from the accumulator, and they
differ by about a factor sqrt(N):

    --noise-map / --variance-map           noise (or variance) of the STACKED
        MEAN, i.e. the uncertainty left in the exported master. This is the
        one to use for the master's SNR. Pixels seen by a single frame have
        no measurable scatter and are written as NaN; pixels with no data are 0.
    --sample-noise-map / --sample-variance-map   frame-to-frame scatter of ONE
        sub-frame about the mean (what --noise-map used to export before it
        was corrected; kept for variability studies).

All of them come from the profile's Section 10 statistics, with W = LSF_WEIGHT,
M2 = LSF_M2 and W2 = LSF_W2 at each pixel:

    n_eff = W^2 / W2            s^2 = M2 / (W - W2/W)       var(mean) = s^2 * W2 / W^2

so n_eff is exact (it no longer depends on the count plane, which saturates at
65535, nor on an approximated weighting efficiency). --neff-map exports n_eff
itself: an effective-integration map, in frames. A file written before LSF_W2
existed can still export its master and count map, but not these maps.

USAGE:
    pyscript lsf_export.py --input stack1.lsf.fits --output stack1.fit
    pyscript lsf_export.py --input stack1.lsf.fits --output stack1_noise.fit --noise-map
    pyscript lsf_export.py --input stack1.lsf.fits --output stack1_var.fit --variance-map
    pyscript lsf_export.py --input stack1.lsf.fits --output stack1_scatter.fit --sample-noise-map
    pyscript lsf_export.py --input stack1.lsf.fits --output stack1_neff.fit --neff-map

Run with no arguments from the Scripts menu (GUI mode) for a small Tkinter
file-picker form instead of the command line.
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lsf_common as lc  # noqa: E402

import numpy as np  # noqa: E402

try:
    import sirilpy as s
    from sirilpy import LogColor
except ImportError:
    print("This script must be run from within Siril (sirilpy not found).")
    raise


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="lsf_export.py",
        description="Export an LSF file's flattened master (or noise/variance map) to plain FITS.",
    )
    p.add_argument("--input", dest="input", default=None, required=False,
                    help="Path to the LSF FITS file to export from.")
    p.add_argument("--output", dest="output", default=None, required=False,
                    help="Output file path.")
    group = p.add_mutually_exclusive_group()
    group.add_argument("--noise-map", dest="noise_map", action="store_true", default=False,
                        help="Export the per-pixel noise (standard error) of the stacked mean: "
                             "sqrt(s^2 * W2 / W^2), s^2 = M2 / (W - W2/W). Single-frame pixels are NaN.")
    group.add_argument("--variance-map", dest="variance_map", action="store_true", default=False,
                        help="Export the per-pixel variance of the stacked mean (the square of --noise-map).")
    group.add_argument("--sample-noise-map", dest="sample_noise_map", action="store_true", default=False,
                        help="Export sqrt(s^2): the unbiased scatter of ONE frame about the mean "
                             "(this was --noise-map before it was corrected; about sqrt(N) larger "
                             "than the noise of the stacked mean). Single-frame pixels are NaN.")
    group.add_argument("--sample-variance-map", dest="sample_variance_map", action="store_true", default=False,
                        help="Export s^2 = M2 / (W - W2/W): the unbiased per-frame variance about "
                             "the mean. Single-frame pixels are NaN.")
    group.add_argument("--neff-map", dest="neff_map", action="store_true", default=False,
                        help="Export the effective sample size n_eff = W^2 / W2 per pixel "
                             "(an effective-integration map, in frames).")
    group.add_argument("--count-map", dest="count_map", action="store_true", default=False,
                        help="Export the per-pixel contribution count instead of the master.")
    return p


def run_gui_and_get_args() -> "argparse.Namespace | None":
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox

    result = {}
    root = tk.Tk()
    root.title("LSF Export")

    input_var = tk.StringVar(value="")
    output_var = tk.StringVar(value="")
    what_var = tk.StringVar(value="master")

    def pick_input():
        path = filedialog.askopenfilename(title="LSF file to export",
                                           filetypes=[("LSF FITS", "*.fits *.fit"), ("All files", "*.*")])
        if path:
            input_var.set(path)

    def pick_output():
        path = filedialog.asksaveasfilename(title="Output file", defaultextension=".fit")
        if path:
            output_var.set(path)

    row = 0
    ttk.Label(root, text="LSF input file:").grid(row=row, column=0, sticky="w", padx=4, pady=2)
    in_frame = ttk.Frame(root)
    ttk.Entry(in_frame, textvariable=input_var, width=30).pack(side="left")
    ttk.Button(in_frame, text="...", width=3, command=pick_input).pack(side="left")
    in_frame.grid(row=row, column=1, sticky="we", padx=4, pady=2)
    row += 1

    ttk.Label(root, text="Export:").grid(row=row, column=0, sticky="w", padx=4, pady=2)
    ttk.Combobox(root, textvariable=what_var, values=["master", "noise_map", "variance_map", "sample_noise_map",
                         "sample_variance_map", "neff_map", "count_map"],
                 state="readonly", width=27).grid(row=row, column=1, sticky="we", padx=4, pady=2)
    row += 1

    ttk.Label(root, text="Output file:").grid(row=row, column=0, sticky="w", padx=4, pady=2)
    out_frame = ttk.Frame(root)
    ttk.Entry(out_frame, textvariable=output_var, width=30).pack(side="left")
    ttk.Button(out_frame, text="...", width=3, command=pick_output).pack(side="left")
    out_frame.grid(row=row, column=1, sticky="we", padx=4, pady=2)
    row += 1

    def on_export():
        if not input_var.get().strip() or not output_var.get().strip():
            messagebox.showerror("LSF Export", "Please choose both an input and output file.")
            return
        result["input"] = input_var.get().strip()
        result["output"] = output_var.get().strip()
        result["noise_map"] = what_var.get() == "noise_map"
        result["variance_map"] = what_var.get() == "variance_map"
        result["sample_noise_map"] = what_var.get() == "sample_noise_map"
        result["sample_variance_map"] = what_var.get() == "sample_variance_map"
        result["neff_map"] = what_var.get() == "neff_map"
        result["count_map"] = what_var.get() == "count_map"
        root.destroy()

    def on_cancel():
        root.destroy()

    btns = ttk.Frame(root)
    ttk.Button(btns, text="Export", command=on_export).pack(side="left", padx=4)
    ttk.Button(btns, text="Cancel", command=on_cancel).pack(side="left", padx=4)
    btns.grid(row=row, column=0, columnspan=2, pady=8)

    root.mainloop()

    if not result:
        return None
    return argparse.Namespace(input=result["input"], output=result["output"],
                                noise_map=result["noise_map"], variance_map=result["variance_map"],
                                sample_noise_map=result["sample_noise_map"],
                                sample_variance_map=result["sample_variance_map"],
                                neff_map=result["neff_map"], count_map=result["count_map"])


# Selectable error maps: flag -> (lsf_common.error_map kind, is a variance, FITS comment, log label).
MAP_EXPORTS = {
    "noise_map": ("mean_noise", False, "std error of the stacked mean per pixel",
                  "noise map (standard error of the stacked mean)"),
    "variance_map": ("mean_variance", True, "variance of the stacked mean per pixel",
                     "variance map (variance of the stacked mean)"),
    "sample_noise_map": ("sample_noise", False, "sqrt(s^2): single-frame scatter per pixel",
                         "sample noise map (single-frame scatter)"),
    "sample_variance_map": ("sample_variance", True, "s^2=M2/(W-W2/W): single-frame variance per pixel",
                            "sample variance map (single-frame variance)"),
}


# Keywords that identify a file as LSF (FITS profile Section 4.2).
LSF_IDENTITY_KEYS = ("LSF", "LSFVER", "LSFUID", "LSFGEN", "LSFSENS", "LSFBAYR", "LSFDMOD", "LSFNORM",
                     "LSFNMIN", "LSFNMAX", "LSFDSUM", "LSFBUSY", "LSFCSAT")


def _plain_keywords(keywords: list, data: dict) -> list:
    """Stack header cards for an exported plain image: identity keywords removed, provenance added.
    The inverse of output normalization is kept as HISTORY, because LSFNMIN / LSFNMAX go with the
    identity keywords."""
    drop = set(LSF_IDENTITY_KEYS)
    kws = [k for k in keywords if k[0].upper() not in drop]
    meta = data.get("meta") or {}
    if meta.get("uid"):
        kws.append(("HISTORY", f"exported from LSF file {meta['uid']}, generation {meta.get('gen')}", ""))
    # "Copy HDU 0": ROWORDER is part of how the pixels must be read (Siril assumes BOTTOM-UP when it
    # is absent, which would flip a TOP-DOWN stack), so it must survive into the plain file.
    if meta.get("row_order"):
        kws.append(("ROWORDER", meta["row_order"], "row order of the pixel array"))
    norm = lc.loads(data["properties"]["LSF:OutputNorm"])
    if norm.get("enabled"):
        kws.append(("HISTORY", "Stack was output-normalized: raw = value * (max - min) + min, "
                    f"min={norm['mini']!r} max={norm['maxi']!r}", ""))
    return kws


def run(args: argparse.Namespace, siril: "s.SirilInterface | None") -> None:
    data = lc.read_lsf(args.input)
    images = data["images"]

    selected = [k for k in list(MAP_EXPORTS) + ["neff_map", "count_map"] if getattr(args, k, False)]
    sel = selected[0] if selected else None
    # The master and the count map do not need LSF_W2, so a file from before that plane existed can
    # still export them; every noise / variance / n_eff map does. Not a validation (C1).
    lc.check_lsf_file(data, args.input, optional_images=() if sel not in (None, "count_map") else ("SumW2",))

    stack = images["Stack"]["data"]
    weight = images["WeightSum"]["data"]
    m2 = images["SumSq"]["data"]
    w2 = images["SumW2"]["data"] if "SumW2" in images else None
    count = images["ContributionCount"]["data"]
    fits_keywords = _plain_keywords(images["Stack"]["fits_keywords"], data)
    # Noise/variance are stored on the raw scale; if Stack is output-normalized, express the
    # maps on the same scale as the exported master (divide by R, or R^2 for variance).
    _norm = lc.loads(data["properties"]["LSF:OutputNorm"])
    norm_range = (float(_norm["maxi"]) - float(_norm["mini"])) if _norm.get("enabled") else 1.0

    if sel in MAP_EXPORTS:
        kind, is_variance, kw_comment, label = MAP_EXPORTS[sel]
        out_data = lc.error_map(kind, m2, weight, w2)
        # Stack may be output-normalized; put the map on the same scale as the exported master.
        out_data /= np.float32(norm_range ** 2 if is_variance else norm_range)
        fits_keywords.append(("LSFEXPRT", sel, kw_comment))
        fits_keywords.append(("HISTORY", "derived from LSF_WEIGHT (W), LSF_M2 (M2), LSF_W2 (W2): "
                              "n_eff = W^2/W2, s^2 = M2/(W - W2/W)"
                              + (", var(mean) = s^2*W2/W^2" if kind.startswith("mean") else "")
                              + "; NaN where n_eff <= 1, 0 where no data", ""))
    elif sel == "neff_map":
        out_data = lc.error_map("n_eff", m2, weight, w2)       # frames, not flux: no output-norm rescale
        fits_keywords.append(("LSFEXPRT", "neff_map", "effective sample size W^2/W2 per pixel"))
        fits_keywords.append(("HISTORY", "n_eff = LSF_WEIGHT^2 / LSF_W2; 0 where no data", ""))
        label = "effective-sample-size (n_eff) map"
    elif sel == "count_map":
        out_data = np.asarray(count, dtype=np.float32)
        fits_keywords.append(("LSFEXPRT", "count_map", "per-pixel contribution count"))
        label = "contribution-count map"
    else:
        out_data = np.asarray(stack, dtype=np.float32)
        label = "flattened master"

    try:
        lc.write_fits(args.output, out_data, fits_keywords)
        if siril is not None:
            siril.log("FITS written directly (astropy)", LogColor.GREEN)
    except (ImportError, RuntimeError) as e:
        if siril is None:
            raise RuntimeError("FITS export requires either astropy or a connected Siril instance.")
        siril.log(f"Falling back to sirilpy for FITS export (astropy unavailable: {e})", LogColor.SALMON)
        header_str = lc.fits_header_string(fits_keywords)
        try:
            ok = siril.save_image_file(out_data, header_str, args.output, timeout=300000)
        except TypeError:
            ok = siril.save_image_file(out_data, header_str, args.output)
        if not ok:
            raise RuntimeError("Siril reported save_image_file() failed.")

    # Release the memory maps of the input file.
    del stack, weight, m2, w2, count, images, data

    msg = f"Exported {label} from '{args.input}' to '{args.output}'."
    if siril is not None:
        siril.log(msg, LogColor.GREEN)
    else:
        print(msg)


def main():
    siril = s.SirilInterface()
    siril.connect()
    try:
        parser = build_arg_parser()
        argv = sys.argv[1:]
        if not argv and not siril.is_cli():
            args = run_gui_and_get_args()
            if args is None:
                siril.log("LSF export cancelled by user.")
                return
        else:
            args = parser.parse_args(argv)
            if not args.input or not args.output:
                raise RuntimeError("--input and --output are required (or use the GUI form).")
        run(args, siril)
    except Exception as e:
        siril.log(f"LSF export failed: {e}", LogColor.RED)
        traceback.print_exc()
        try:
            siril.error_messagebox(f"LSF export failed:\n{e}")
        except Exception:
            pass
        sys.exit(1)
    finally:
        siril.disconnect()


if __name__ == "__main__":
    main()