#!/usr/bin/env python3
"""
lsf_stack.py -- Living Stack Format (LSF) incremental stacker for Siril.

Author: (fill in your name here before distributing)
License: (fill in a license here before distributing, e.g. GPL-3.0)
Contact: (fill in a contact / forum / repo URL here before distributing)

WHAT THIS SCRIPT DOES
----------------------
Takes a Siril sequence whose frames are already calibrated, registered and
plate-solved (this script performs NO calibration, registration or plate
solving of its own -- see live_stack.md, Section 4 "Explicitly Out of
Scope" and the "final stacking only" brief), and folds it into a Living
Stack Format accumulator: a compact set of per-pixel sufficient statistics
(weighted mean, Welford M2, sum of weights, sum of squared weights, contribution
count) plus a frame ledger,
stored in an LSF file (a FITS file, see live_stack_fits.md).

You can either:
  - start a brand new LSF file from this sequence, or
  - fold this sequence into an *existing* LSF file's history (--append),
    or
  - stack this sequence fresh, ignoring any existing LSF file's history
    even if one is given (default when --append is not passed).

Rejection and weighting are modelled on Siril's own "stack ... rej ..."
command -- see --rejection / --weighting below and lsf_common.py for
exactly which of Siril's algorithms are reproduced and which are
deliberately out of scope for a bounded accumulator.

Normalization is NOT a menu of interchangeable options: live_stack.md
Section 8.3 mandates one specific ingestion-time algorithm (a robust
linear fit of each incoming frame against the accumulator's own current
best estimate, or against the best-quality frame in the batch when there
is no existing accumulator yet), and this script always applies it unless
--skip-normalization is passed. See lsf_common.py's fit_normalization().

Mono and one-shot-color (OSC) sensors are both supported (Section 10):
normalization, weight, and background are tracked independently per
color plane. A file's sensor class (mono / CFA + Bayer pattern), filter,
and demosaic method are fixed at creation and enforced on every later
--append (Section 10.1 / 9).

USAGE (from Siril's command line / a script, or the Scripts menu):
    pyscript lsf_stack.py --seq lightseq --output stack1.lsf.fits
    pyscript lsf_stack.py --seq lightseq --existing stack1.lsf.fits \\
        --append --output stack1.lsf.fits \\
        --rejection winsorized --sigma-low 3 --sigma-high 3 \\
        --weighting noise --output-norm

Run with no arguments from the Scripts menu (GUI mode) to get a small
Tkinter form instead of the command line.

CONTAINER: the output is an LSF FITS file (live_stack_fits.md, v0.3). On every append the existing
file is validated first (Section 8) and a failed check refuses the merge and leaves the file
untouched. The new file is written beside the old one and renamed over it (Section 10, preferred
commit protocol).
"""

from __future__ import annotations

import argparse
import gc
import os
import re
import shutil
import sys
import tempfile
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lsf_common as lc  # noqa: E402

import numpy as np  # noqa: E402

try:
    import sirilpy as s
    from sirilpy import LogColor
except ImportError:
    print("This script must be run from within Siril (sirilpy not found).")
    raise


# FITS profile Section 10 (Initialize): warn when the first batch has fewer frames than this.
LOW_STACK_WARN_FRAMES = 25

SUPPORTED_SENSOR_CLASS = ("mono", "cfa")
SUPPORTED_BAYER_PATTERN = ("RGGB", "GRBG", "GBRG", "BGGR", "none", "unknown")
SUPPORTED_DEMOSAIC_METHOD = ("mono", "debayer_first", "cfa_drizzle")


# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="lsf_stack.py",
        description="Stack a Siril sequence into a Living Stack Format (LSF) file.",
    )
    p.add_argument("--seq", dest="seq", default=None,
                    help="Sequence name to stack. If omitted, the currently loaded sequence is used.")
    p.add_argument("--existing", dest="existing", default=None,
                    help="Path to an existing LSF FITS file whose history this batch may extend.")
    p.add_argument("--append", dest="append", action="store_true", default=False,
                    help="Fold this batch into --existing's accumulated history. "
                         "Without this flag, --existing (if given) is ignored and a fresh, "
                         "independent stack is produced from this batch alone.")
    p.add_argument("--output", dest="output", default=None,
                    help="Output LSF path (.lsf.fits). May be the same "
                         "path as --existing for true in-place growth, or a different path to keep a "
                         "versioned snapshot.")
    p.add_argument("--keep-backup", dest="keep_backup", action="store_true", default=False,
                    help="When --output already exists (in-place append), keep the previous file as "
                         "<output>.bak instead of replacing it (Section 10).")
    p.add_argument("--rejection", dest="rejection", default=None,
                    choices=lc.SUPPORTED_REJECTION,
                    help="Pixel rejection algorithm for this batch (default: sigma, or the existing file's stored value with --append). "
                         "Recommended by frame count: 0-4 frames=percentile, 5-10=sigma, "
                         "11-30=winsorized, 31-300=gesdt, 301+=linearfit. "
                         "Note: these work within the current batch only, not across accumulated history.")
    p.add_argument("--sigma-low", dest="sigma_low", type=float, default=None,
                    help="Low threshold; default is per-method (sigma 3, linearfit 5, percentile 0.2).")
    p.add_argument("--sigma-high", dest="sigma_high", type=float, default=None,
                    help="High threshold; default is per-method (sigma 3, linearfit 5, percentile 0.1).")
    p.add_argument("--max-iters", dest="max_iters", type=int, default=None,
                    help="Maximum number of rejection passes within this batch (default: 3). "
                         "Each pass can reject up to (samples - 4) per pixel; a pixel never "
                         "drops below 4 samples.")
    p.add_argument("--gesdt-alpha", dest="gesdt_alpha", type=float, default=None,
                    help="GESDT significance level (default: 0.05 for 95%% confidence). "
                         "Used only when --rejection=gesdt.")
    p.add_argument("--gesdt-max-frac", dest="gesdt_max_frac", type=float, default=None,
                    help="GESDT maximum outlier fraction per pixel (default: 0.3).")
    p.add_argument("--weighting", dest="weighting", default=None,
                    choices=lc.SUPPORTED_WEIGHTING)
    p.add_argument("--skip-normalization", dest="skip_normalization", action=argparse.BooleanOptionalAction, default=None,
                    help="Debug/testing escape hatch: use scale=1, offset=0 for every frame instead "
                         "of the Section 8.3 linear-fit normalization. Not the spec's default "
                         "behavior -- only use this if you specifically need to bypass it.")
    p.add_argument("--saturation-level", dest="saturation_level", type=float, default=None,
                    help="Exclude pixels at/above this value from the normalization fit's overlap "
                         "region (bit-depth-agnostic; off by default). Does not affect stacking itself.")
    p.add_argument("--sensor-class", dest="sensor_class", default=None, choices=SUPPORTED_SENSOR_CLASS,
                    help="Declared sensor class (Section 10.1). Auto-detected from channel count "
                         "if omitted (1 -> mono, 3 -> cfa). Enforced on every later --append.")
    p.add_argument("--bayer-pattern", dest="bayer_pattern", default=None, choices=SUPPORTED_BAYER_PATTERN,
                    help="Bayer pattern, only meaningful when --sensor-class=cfa. Defaults to "
                         "'unknown' (recorded, but not validated as strictly) if not given.")
    p.add_argument("--demosaic-method", dest="demosaic_method", default=None, choices=SUPPORTED_DEMOSAIC_METHOD,
                    help="Fixed at file creation (Section 10.3). Auto: 'mono' for mono sensors, "
                         "'debayer_first' for CFA sensors if not given.")
    p.add_argument("--norm-method", dest="norm_method", default=None, choices=list(lc.NORM_METHODS),
                    help="Ingestion-time normalization (fixed at file creation). linearfit: robust "
                         "signal fit against the accumulator mean (Section 8.3, default). addscale: "
                         "Siril's -norm=addscale, matching each frame's background location and "
                         "spread to a stored baseline frame (parity mode).")
    p.add_argument("--norm-estimator", dest="norm_estimator", default=None, choices=list(lc.NORM_ESTIMATORS),
                    help="addscale only. lite: median and 1.5*MAD computed here on the merged pixels "
                         "(Siril -lite_norm). siril: Siril's IKSS location/scale from the frame statistics "
                         "(measured before alignment). Default: lite.")
    p.add_argument("--equalize-rgb", dest="equalize_rgb", action=argparse.BooleanOptionalAction, default=None,
                    help="Color files: match every channel's background location and spread to the green "
                         "channel of the reference frame (Siril -rgb_equal; this is what makes the "
                         "per-channel noise of an OSC stack roughly equal). Works with both "
                         "normalization methods. Fixed at file creation. Default: on for 3-channel "
                         "files, off for mono; --no-equalize-rgb disables it.")
    p.add_argument("--ref-frame", dest="ref_frame", type=int, default=None,
                    help="On initialization: 1-based sequence frame number to use as the normalization "
                         "reference / baseline (Siril's registration reference). Default: the "
                         "best-quality frame. Ignored on append (the stored baseline is used).")
    p.add_argument("--output-norm", dest="output_norm", action=argparse.BooleanOptionalAction, default=None,
                    help="Store the Stack plane min-max normalized to [0,1] (Siril norm_to_0_1_range). "
                         "The transform is recorded in the file and inverted on --append, so the "
                         "accumulator stays on the raw scale. Use --no-output-norm to store the "
                         "raw mean. Default: inherited from --existing with --append, else off.")
    p.add_argument("--skip-align", dest="skip_align", action="store_true", default=False,
                    help="With --append: do NOT align the new frames to the existing LSF framing; "
                         "just add them to the accumulator as they are (they must already share "
                         "its pixel grid). Default: align via Siril registration.")
    p.add_argument("--align-transform", dest="align_transform", default=None,
                    choices=["shift", "similarity", "affine", "homography"],
                    help="Transformation Siril registers new frames with (default: homography, or "
                         "the existing file's stored value).")
    p.add_argument("--align-interp", dest="align_interp", default=None,
                    choices=["nearest", "linear", "cubic", "lanczos4", "area"],
                    help="Interpolation used when resampling onto the LSF framing (default: lanczos4).")
    p.add_argument("--align-minpairs", dest="align_minpairs", type=int, default=None,
                    help="Minimum matched star pairs for a frame to be accepted; frames below it "
                         "are dropped, not stacked (default: 10).")
    p.add_argument("--compress", dest="compress", action="store_true", default=False,
                    help="Tile-compress the extension images (LSF_WEIGHT / LSF_M2 / LSF_W2: GZIP_2, lossless; "
                         "LSF_COUNT: RICE_1). Off by default for maximum compatibility. The primary "
                         "Stack is never compressed (C1). Expect roughly 20-35 percent smaller files "
                         "and slower appends.")
    p.add_argument("--min-frames-warn", dest="min_frames_warn", type=int, default=LOW_STACK_WARN_FRAMES,
                    help="When a NEW LSF file is created from fewer than this many frames, warn before "
                         "writing it (FITS profile Section 10): per-pixel scatter estimated from so few "
                         "frames is unreliable, and any rejection decision made on it stays in the "
                         "stack. The warning never blocks creation. 0 disables it "
                         f"(default: {LOW_STACK_WARN_FRAMES}).")
    p.add_argument("--prefetch", dest="prefetch", default="spill", choices=["off", "spill", "all"],
                    help="Overlap reads with compute using one look-ahead worker thread. spill "
                         "(default): Pass 2 reads the next strip of the normalized-frame spill and of "
                         "the existing accumulator while the current strip is being rejected and "
                         "combined. all: also reads the next Siril frame during Pass 1 (sirilpy "
                         "serializes commands with a lock, but this has not been run against a real "
                         "Siril; use off or spill if a script run misbehaves). off: no threads. "
                         "Results are identical in every mode; the look-ahead costs about one more "
                         "strip of RAM, which is accounted for in the strip size.")
    p.add_argument("--mem-budget-gb", dest="mem_budget_gb", type=float, default=2.0,
                    help="Working-memory budget (GB) for the rejection/combine strips "
                         "(default: 2.0). Peak RAM is roughly this plus the accumulator planes.")
    p.add_argument("--temp-dir", dest="temp_dir", default=None,
                    help="Folder for the temporary normalized-frame spill file (about N x one "
                         "float32 frame of disk). Default: an 'lsf_temp' folder next to --output. "
                         "The spill folder is always deleted when processing finishes or fails.")
    p.add_argument("--device", dest="device", default="auto",
                    choices=["auto", "cpu", "cuda", "mps", "off"],
                    help="Compute backend for rejection/combine. auto = PyTorch on CUDA/MPS, else "
                         "PyTorch CPU tensors; off = NumPy only. On a PyTorch out-of-memory error the "
                         "strip is split and retried; any other PyTorch error falls back to NumPy "
                         "for that strip, with a warning (default: auto).")
    p.add_argument("--torch-precision", dest="torch_precision", default="auto",
                    choices=["auto", "float32", "float64"],
                    help="dtype of the PyTorch rejection work (sort, median, thresholds). auto = "
                         "float32 on CUDA and MPS, float64 on CPU tensors. The final weighted sums "
                         "stay float64 wherever the device supports it. float64 on a consumer GPU "
                         "is 30-60x slower; use it on datacenter GPUs if you want bit-for-bit "
                         "agreement with the NumPy path.")
    p.add_argument("--object-name", dest="object_name", default=None,
                    help="Optional OBJECT name to embed in the output header.")
    return p


# --------------------------------------------------------------------------
# Small Tkinter fallback GUI, used only when the script is launched with no
# arguments from Siril's Scripts menu (i.e. not from the command line).
# --------------------------------------------------------------------------

def run_gui_and_get_args(default_seqname: str) -> "argparse.Namespace | None":
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox

    result = {}
    root = tk.Tk()
    root.title("LSF Stack")

    row = 0

    def add_row(label, widget_factory):
        nonlocal row
        ttk.Label(root, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=2)
        widget = widget_factory()
        widget.grid(row=row, column=1, sticky="we", padx=4, pady=2)
        row += 1
        return widget

    seq_var = tk.StringVar(value=default_seqname or "")
    add_row("Sequence name:", lambda: ttk.Entry(root, textvariable=seq_var, width=32))

    existing_var = tk.StringVar(value="")

    def pick_existing():
        path = filedialog.askopenfilename(title="Existing LSF file (optional)",
                                           filetypes=[("LSF FITS", "*.fits *.fit"), ("All files", "*.*")])
        if path:
            existing_var.set(path)
            try:
                info = lc.read_lsf_settings(path)
            except Exception as e:
                messagebox.showwarning("LSF Stack", f"Could not read settings from file:\n{e}")
                return
            st = info["settings"]
            append_var.set(True)
            rej_var.set(st.get("rejection") or "sigma")
            sigl_var.set("" if st.get("sigma_low") is None else str(st["sigma_low"]))
            sigh_var.set("" if st.get("sigma_high") is None else str(st["sigma_high"]))
            weight_var.set(st.get("weighting") or "none")
            outnorm_var.set(bool(st.get("output_norm")))
            skipnorm_var.set(bool(st.get("skip_normalization")))
            normm_var.set(st.get("norm_method") or "linearfit")
            esti_var.set(st.get("norm_estimator") or "lite")
            eqrgb_var.set(bool(st.get("equalize_rgb")))
            messagebox.showinfo("LSF Stack", f"Loaded stored settings from the file "
                                f"({info['n_frames']} frames in history).")

    existing_frame = ttk.Frame(root)
    ttk.Entry(existing_frame, textvariable=existing_var, width=26).pack(side="left")
    ttk.Button(existing_frame, text="...", width=3, command=pick_existing).pack(side="left")
    ttk.Label(root, text="Existing LSF (optional):").grid(row=row, column=0, sticky="w", padx=4, pady=2)
    existing_frame.grid(row=row, column=1, sticky="we", padx=4, pady=2)
    row += 1

    append_var = tk.BooleanVar(value=False)
    ttk.Checkbutton(root, text="Append to existing (else: fresh, independent stack)",
                     variable=append_var).grid(row=row, column=0, columnspan=2, sticky="w", padx=4)
    row += 1

    output_var = tk.StringVar(value="")

    def pick_output():
        path = filedialog.asksaveasfilename(title="Output LSF file", defaultextension=".lsf.fits",
                                             filetypes=[("LSF FITS", "*.fits")])
        if path:
            output_var.set(path)

    output_frame = ttk.Frame(root)
    ttk.Entry(output_frame, textvariable=output_var, width=26).pack(side="left")
    ttk.Button(output_frame, text="...", width=3, command=pick_output).pack(side="left")
    ttk.Label(root, text="Output LSF file:").grid(row=row, column=0, sticky="w", padx=4, pady=2)
    output_frame.grid(row=row, column=1, sticky="we", padx=4, pady=2)
    row += 1

    rej_var = tk.StringVar(value="sigma")
    add_row("Rejection:", lambda: ttk.Combobox(root, textvariable=rej_var,
                                                values=list(lc.SUPPORTED_REJECTION), state="readonly", width=29))
    sigl_var = tk.StringVar(value="")
    add_row("Sigma low (blank=default):", lambda: ttk.Entry(root, textvariable=sigl_var, width=32))
    sigh_var = tk.StringVar(value="")
    add_row("Sigma high (blank=default):", lambda: ttk.Entry(root, textvariable=sigh_var, width=32))
    weight_var = tk.StringVar(value="none")
    add_row("Weighting:", lambda: ttk.Combobox(root, textvariable=weight_var,
                                                 values=list(lc.SUPPORTED_WEIGHTING), state="readonly", width=29))
    sensor_var = tk.StringVar(value="(auto)")
    add_row("Sensor class:", lambda: ttk.Combobox(root, textvariable=sensor_var,
                                                    values=["(auto)"] + list(SUPPORTED_SENSOR_CLASS),
                                                    state="readonly", width=29))
    normm_var = tk.StringVar(value="linearfit")
    add_row("Normalization method:", lambda: ttk.Combobox(root, textvariable=normm_var,
                                                           values=list(lc.NORM_METHODS), state="readonly", width=29))
    esti_var = tk.StringVar(value="lite")
    add_row("addscale estimator:", lambda: ttk.Combobox(root, textvariable=esti_var,
                                                         values=list(lc.NORM_ESTIMATORS), state="readonly", width=29))
    ref_var = tk.StringVar(value="")
    add_row("Reference frame # (blank=auto):", lambda: ttk.Entry(root, textvariable=ref_var, width=32))
    eqrgb_var = tk.BooleanVar(value=True)
    ttk.Checkbutton(root, text="Equalize RGB to green (color only; fixed at creation)",
                     variable=eqrgb_var).grid(row=row, column=0, columnspan=2, sticky="w", padx=4)
    row += 1
    outnorm_var = tk.BooleanVar(value=False)
    ttk.Checkbutton(root, text="Output normalization", variable=outnorm_var).grid(
        row=row, column=0, columnspan=2, sticky="w", padx=4)
    row += 1
    skipnorm_var = tk.BooleanVar(value=False)
    ttk.Checkbutton(root, text="Skip normalization (debug only -- not the spec default)",
                     variable=skipnorm_var).grid(row=row, column=0, columnspan=2, sticky="w", padx=4)
    row += 1
    skipalign_var = tk.BooleanVar(value=False)
    ttk.Checkbutton(root, text="Skip alignment to existing LSF (add frames as-is)",
                     variable=skipalign_var).grid(row=row, column=0, columnspan=2, sticky="w", padx=4)
    row += 1

    def on_run():
        if not seq_var.get().strip():
            messagebox.showerror("LSF Stack", "Please enter a sequence name.")
            return
        if not output_var.get().strip():
            messagebox.showerror("LSF Stack", "Please choose an output file.")
            return
        result["seq"] = seq_var.get().strip()
        result["existing"] = existing_var.get().strip() or None
        result["append"] = bool(append_var.get())
        result["output"] = output_var.get().strip()
        result["rejection"] = rej_var.get()
        result["sigma_low"] = _opt_float(sigl_var.get())
        result["sigma_high"] = _opt_float(sigh_var.get())
        result["weighting"] = weight_var.get()
        result["sensor_class"] = None if sensor_var.get() == "(auto)" else sensor_var.get()
        result["output_norm"] = bool(outnorm_var.get())
        inherit = bool(append_var.get()) and bool(existing_var.get().strip())   # appending: file's values win
        result["norm_method"] = None if inherit else normm_var.get()
        result["norm_estimator"] = None if inherit else esti_var.get()
        result["equalize_rgb"] = None if inherit else bool(eqrgb_var.get())
        try:
            result["ref_frame"] = int(ref_var.get()) if ref_var.get().strip() else None
        except ValueError:
            messagebox.showerror("LSF Stack", "Reference frame must be a whole number.")
            result.clear()
            return
        result["skip_normalization"] = bool(skipnorm_var.get())
        result["skip_align"] = bool(skipalign_var.get())
        root.destroy()

    def on_cancel():
        root.destroy()

    btns = ttk.Frame(root)
    ttk.Button(btns, text="Run", command=on_run).pack(side="left", padx=4)
    ttk.Button(btns, text="Cancel", command=on_cancel).pack(side="left", padx=4)
    btns.grid(row=row, column=0, columnspan=2, pady=8)

    root.mainloop()

    if not result:
        return None
    ns = argparse.Namespace(
        seq=result["seq"], existing=result["existing"], append=result["append"],
        output=result["output"], rejection=result["rejection"], sigma_low=result["sigma_low"],
        sigma_high=result["sigma_high"], max_iters=None, gesdt_alpha=None, gesdt_max_frac=None,
        weighting=result["weighting"],
        skip_normalization=result["skip_normalization"], saturation_level=None,
        skip_align=result["skip_align"], align_transform=None, align_interp=None, align_minpairs=None,
        sensor_class=result["sensor_class"], bayer_pattern=None, demosaic_method=None,
        output_norm=result["output_norm"], compress=False, object_name=None, keep_backup=False,
        mem_budget_gb=2.0, temp_dir=None, device="auto", torch_precision="auto", prefetch="spill",
        min_frames_warn=LOW_STACK_WARN_FRAMES,
        norm_method=result["norm_method"], norm_estimator=result["norm_estimator"],
        equalize_rgb=result["equalize_rgb"], ref_frame=result["ref_frame"],
    )
    return ns


# --------------------------------------------------------------------------
# Siril data extraction helpers
# --------------------------------------------------------------------------

def _to_hwc(arr: np.ndarray) -> np.ndarray:
    """Normalize a pixel array from Siril's (channels, height, width)
    convention (per get_image_shape()'s documented tuple order) to this
    script's internal channel-last (height, width, channels) convention.
    Mono arrays (2D) pass through unchanged, and a single-plane 3D array,
    (1, H, W) or (H, W, 1), is squeezed to (H, W): the rest of the script treats
    mono as 2D, and a trailing length-1 axis breaks the block-median fit.

    If your Siril/sirilpy version actually returns channel-last arrays
    already, this function will mis-transpose color data -- check a known
    RGB frame's shape once (compare against Sequence.rx/ry/nb_layers) and
    flip this if needed; everything downstream only assumes channel-last.
    """
    if arr.ndim == 2:
        return arr
    if arr.ndim == 3:
        c, h, w = arr.shape
        # Single plane: (1, H, W) channel-first or (H, W, 1) channel-last -> (H, W).
        if c == 1 and w != 1:
            return arr[0]
        if w == 1 and c not in (1, 3):
            return arr[..., 0]
        # Heuristic: channel count for astronomical sequences is 1 or 3;
        # if the first axis isn't a plausible channel count but the last
        # one is, assume the array is already channel-last.
        if c == 3 and w not in (1, 3):
            return np.transpose(arr, (1, 2, 0))
        return arr
    raise ValueError(f"unexpected pixel array ndim={arr.ndim}")


def _extract_wcs(header: dict) -> dict:
    keys = ["CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
            "CDELT1", "CDELT2", "CROTA2", "CTYPE1", "CTYPE2", "EQUINOX", "PLTSOLVD"]
    return {k: header[k] for k in keys if k in header}


def _wcs_pixel_scale(wcs: dict) -> "float | None":
    """Rough pixel scale in degrees/pixel, for a sanity check that two
    grids are compatible -- not used for anything else."""
    if "CD1_1" in wcs and "CD2_1" in wcs:
        return float(np.hypot(wcs["CD1_1"], wcs["CD2_1"]))
    if "CDELT1" in wcs:
        return abs(float(wcs["CDELT1"]))
    return None


def _extract_filter(header: dict) -> "str | None":
    val = header.get("FILTER")
    if val is None:
        return None
    val = str(val).strip()
    return val or None


def _infer_sensor_profile(channels: int, header: dict, args: argparse.Namespace) -> dict:
    """
    Resolve this batch's sensor_class / bayer_pattern / demosaic_method
    (live_stack.md Section 10.1 / 10.3), from explicit CLI args where
    given, else a best-effort auto-detection from channel count and
    whatever Bayer-related FITS keyword the frame happens to carry
    (BAYERPAT is common; not guaranteed present, since by the time this
    script sees the frame it has typically already been debayered or
    Bayer-drizzled by whatever did registration).
    """
    sensor_class = args.sensor_class
    if sensor_class is None:
        sensor_class = "cfa" if channels == 3 else "mono"

    bayer_pattern = args.bayer_pattern
    if bayer_pattern is None:
        if sensor_class == "mono":
            bayer_pattern = "none"
        else:
            raw = str(header.get("BAYERPAT", "")).strip().upper()
            bayer_pattern = raw if raw in SUPPORTED_BAYER_PATTERN else "unknown"

    demosaic_method = args.demosaic_method
    if demosaic_method is None:
        demosaic_method = "mono" if sensor_class == "mono" else "debayer_first"

    return {
        "sensor_class": sensor_class,
        "bayer_pattern": bayer_pattern,
        "demosaic_method": demosaic_method,
    }


def _pick_best_quality_frame(frames: list) -> int:
    """
    Pick the batch-local index of the frame to use as the ingestion-time
    normalization reference when there is no existing accumulator to
    normalize against yet -- live_stack.md Section 8.3: "the first batch
    is normalized ... against the best-quality frame in that batch (chosen
    using the Section 8.1 metrics)". That frame's own photometric scale
    becomes the file's permanent baseline.

    The spec names the available metrics (FWHM, roundness, star count) but
    not a formula for combining them into one ranking, so this is a
    documented heuristic, not a spec-mandated one: more detected stars
    first, tie-broken by tighter FWHM, tie-broken by roundness closest to
    1. Falls back to frame 0 if no frame in the batch has usable
    registration data at all.
    """
    scored = []
    for n, f in enumerate(frames):
        reg = f.get("reg")
        if reg is None:
            continue
        stars = getattr(reg, "number_of_stars", None)
        fwhm = getattr(reg, "fwhm", None)
        roundness = getattr(reg, "roundness", None)
        if stars is None and fwhm is None:
            continue
        scored.append((
            -(stars if stars is not None else 0),
            fwhm if fwhm is not None else float("inf"),
            abs((roundness if roundness is not None else 1.0) - 1.0),
            n,
        ))
    if not scored:
        return 0
    scored.sort()
    return scored[0][-1]


def gather_frame_meta(siril: "s.SirilInterface", idx: int, n_channels: int,
                       channel_for_reg: int = 0) -> dict:
    """Metadata only (header, registration data, per-channel stats) -- NO
    pixel data. Pixels are fetched one frame at a time by read_frame_pixels()
    so memory does not grow with the number of frames."""
    header = siril.get_seq_frame_header(idx, return_as="dict") or {}
    header_str = siril.get_seq_frame_header(idx, return_as="str") or ""
    filename = siril.get_seq_frame_filename(idx)

    reg = None
    try:
        reg = siril.get_seq_regdata(idx, channel_for_reg)
    except Exception:
        reg = None

    # Per-channel stats (Section 10.2: weight/background tracked per plane).
    stats_per_channel = []
    for c in range(n_channels):
        try:
            stats_per_channel.append(siril.get_seq_stats(idx, c))
        except Exception:
            stats_per_channel.append(None)

    return {
        "index": idx,
        "header": header,
        "header_str": header_str,
        "filename": filename,
        "reg": reg,
        "stats_per_channel": stats_per_channel,
    }


def read_frame_pixels(siril: "s.SirilInterface", idx: int, scale: float = 1.0) -> np.ndarray:
    """Read one frame's pixels as channel-last float32 (no copy if already float32)."""
    px = _to_hwc(siril.get_seq_frame_pixeldata(idx, preview=False)).astype(np.float32, copy=False)
    return px * np.float32(scale) if scale != 1.0 else px


# --------------------------------------------------------------------------
# Main pipeline
# --------------------------------------------------------------------------

ALIGN_ROOT = "lsfal"       # temp sequence: lsfal00001.fit (= LSF reference), lsfal00002.fit ...
ALIGN_PREFIX = "a_"        # Siril output prefix -> a_lsfal


def _q(path: str) -> str:
    return f'"{path}"' if " " in path else path


# Per-frame image statistics copied from Siril's per-frame stats into the ledger's
# `quality` block (one list entry per colour plane). bgnoise is stored separately
# as `noise` (and `noise_normalized`); these are the other cheap, useful ones.
LEDGER_STAT_FIELDS = ("mean", "sigma", "mad", "sqrtbwmv", "location", "scale")


def _num(v):
    """float(v), or None when v is missing / not a finite number (keeps the ledger
    JSON clean: numpy scalars would otherwise be written as strings)."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def _h_list(reg):
    H = getattr(reg, "H", None) if reg is not None else None
    if H is None:
        return None
    try:
        return [float(getattr(H, f"h{i}{j}")) for i in range(3) for j in range(3)]
    except Exception:
        return None


def _find_seq(folder: str, name: str):
    for cand in (name, name + "_"):
        if os.path.exists(os.path.join(folder, cand + ".seq")):
            return cand
    return None


class AlignContext:
    """Remembers what to undo after the alignment stage: Siril's working
    directory, the originally loaded sequence, and the temp folder."""

    def __init__(self, tmpdir: str, root: str):
        self.tmpdir, self.root = tmpdir, root
        self.orig_wd = None
        self.orig_seq = None
        self.cd_done = False

    def restore(self, siril, log, success: bool) -> None:
        if self.cd_done and self.orig_wd:
            try:
                siril.cmd("cd", _q(self.orig_wd))
            except Exception as e:
                log(f"Could not restore Siril working directory: {e}", LogColor.SALMON)
        if self.cd_done and self.orig_seq:
            try:
                siril.cmd("load_seq", self.orig_seq)
            except Exception as e:
                log(f"Could not reload original sequence '{self.orig_seq}': {e}", LogColor.SALMON)
        if success:
            shutil.rmtree(self.tmpdir, ignore_errors=True)
            try:
                os.rmdir(self.root)      # only if empty
            except OSError:
                pass
        else:
            log(f"Alignment temp folder kept for inspection: {self.tmpdir}", LogColor.SALMON)


def align_to_lsf(siril, args, frames: list, seqname: str, channels: int, state: dict, log) -> list:
    """Register the new batch onto the existing LSF's framing using Siril.

    1. Export the LSF Stack plane as FITS (frame #1 of a temp sequence, set as the
       registration reference) followed by the selected frames of the input sequence
       (unselected frames are never written, so they are ignored).
    2. `register` the temp sequence to that reference (keeps the LSF framing); fall back
       to `seqapplyreg -framing=current` if register did not write images.
    3. Load the registered sequence and map each surviving frame back to its source frame.
       Frames Siril could not register (below --align-minpairs) are dropped and reported,
       never stacked misaligned.

    Returns the list of surviving frame-meta dicts, re-indexed into the registered sequence,
    with per-frame 'pixel_scale' (65535 if the source was uint16: temp frames are float [0,1]
    and are scaled back on read so the accumulator's units are unchanged).
    """
    ex = lc.read_lsf(args.existing)
    stack = ex["images"]["Stack"]["data"]
    weight = ex["images"]["WeightSum"]["data"]
    ex_ch = stack.shape[2] if stack.ndim == 3 else 1
    if ex_ch != channels:
        raise RuntimeError(f"Cannot align: existing LSF has {ex_ch} channel(s), this batch has {channels}.")
    covered = np.asarray(weight) > 0
    if not covered.any():
        log("Existing LSF accumulator is empty; nothing to align to.")
        del stack, weight, ex, covered
        gc.collect()
        return frames
    ref = np.array(stack, dtype=np.float32)
    ok = covered & np.isfinite(ref)
    lo, hi = float(ref[ok].min()), float(ref[ok].max())
    if hi <= lo:
        hi = lo + 1.0
    ref = np.where(ok, (ref - lo) / (hi - lo), 0.0).astype(np.float32)   # [0,1] for star detection
    if ref.ndim == 3:
        ref = np.ascontiguousarray(np.transpose(ref, (2, 0, 1)))          # Siril channel-first
    del stack, weight, ex, covered, ok
    gc.collect()

    default_tmp = os.path.join(os.path.dirname(os.path.abspath(args.output)), "lsf_temp")
    root = args.temp_dir or default_tmp
    os.makedirs(root, exist_ok=True)
    tmpdir = tempfile.mkdtemp(prefix="lsf_align_", dir=root)
    ctx = AlignContext(tmpdir, root)
    state["ctx"] = ctx
    ctx.orig_wd = siril.get_siril_wd()
    ctx.orig_seq = seqname

    log(f"Aligning {len(frames)} frame(s) to the existing LSF framing "
        f"(transform={args.align_transform}, interp={args.align_interp}, "
        f"min pairs={args.align_minpairs}); temp folder {tmpdir}")

    def save(arr, header, name_idx):
        path = os.path.join(tmpdir, f"{ALIGN_ROOT}{name_idx:05d}.fit")
        if not siril.save_image_file(arr, header, path):
            raise RuntimeError(f"Siril failed to save temporary frame {path}")

    # Step 1: reference first, then the selected frames (read while the original
    # sequence is still the loaded one and the working directory unchanged).
    save(ref, lc.fits_header_string([("OBJECT", "LSF_REFERENCE")]), 1)
    del ref
    scales = []
    for n, f in enumerate(frames):
        px = np.asarray(siril.get_seq_frame_pixeldata(f["index"], preview=False))
        if px.dtype == np.uint16:
            arr, sc = px.astype(np.float32) / 65535.0, 65535.0
        else:
            arr, sc = np.ascontiguousarray(px, dtype=np.float32), 1.0
        scales.append(sc)
        save(arr, f["header_str"] or lc.fits_header_string([("OBJECT", seqname)]), n + 2)
        siril.update_progress(f"Writing alignment frames {n + 1}/{len(frames)}", -1.0)
    n_total = len(frames) + 1

    # Step 2: register to the reference.
    siril.cmd("cd", _q(tmpdir))
    ctx.cd_done = True
    if not siril.create_new_seq(ALIGN_ROOT):
        raise RuntimeError("Siril could not create the temporary alignment sequence.")
    siril.cmd("setref", ALIGN_ROOT, "1")
    siril.update_progress("Registering to LSF reference (Siril)", -1.0)
    siril.cmd("register", ALIGN_ROOT, f"-prefix={ALIGN_PREFIX}", f"-transf={args.align_transform}",
              f"-interp={args.align_interp}", f"-minpairs={args.align_minpairs}")
    out_seq = _find_seq(tmpdir, ALIGN_PREFIX + ALIGN_ROOT)
    if out_seq is None:
        log("register wrote no images; applying the computed registration with seqapplyreg.")
        siril.cmd("seqapplyreg", ALIGN_ROOT, f"-prefix={ALIGN_PREFIX}", f"-interp={args.align_interp}",
                  "-framing=current", "-filter-included")
        out_seq = _find_seq(tmpdir, ALIGN_PREFIX + ALIGN_ROOT)
    if out_seq is None:
        raise RuntimeError("Siril produced no registered sequence (check the Siril log).")

    # Step 3: which frames survived registration, and their transforms.
    siril.cmd("load_seq", ALIGN_ROOT)
    src = siril.get_seq()
    surv = [i for i in range(src.number) if src.imgparam[i].incl]
    if 0 not in surv:
        raise RuntimeError("Siril excluded the LSF reference frame from the temporary sequence.")
    align_h = {}
    for i in surv:
        try:
            align_h[i] = _h_list(siril.get_seq_regdata(i, 0))
        except Exception:
            align_h[i] = None
    dropped = [i for i in range(1, n_total) if i not in surv]

    siril.cmd("load_seq", out_seq)
    aseq = siril.get_seq()
    if aseq.number == len(surv):
        mapping = list(zip(range(aseq.number), surv))      # same order as the surviving frames
    else:                                                  # fall back to the number in the file name
        mapping = []
        for j in range(aseq.number):
            m = re.search(r"(\d+)\.[A-Za-z]+$", os.path.basename(siril.get_seq_frame_filename(j) or ""))
            if not m or not (0 <= int(m.group(1)) - 1 < n_total):
                raise RuntimeError("Cannot map registered frames back to their sources.")
            mapping.append((j, int(m.group(1)) - 1))

    kept = []
    for j, t_idx in mapping:
        if t_idx == 0 or not aseq.imgparam[j].incl:
            continue
        orig = frames[t_idx - 1]
        f = dict(orig)
        f.update(orig_index=orig["index"], orig_filename=orig["filename"], index=j,
                 pixel_scale=scales[t_idx - 1], alignment_H=align_h.get(t_idx),
                 alignment={"method": "siril_register", "transform": args.align_transform,
                            "interp": args.align_interp, "reference": "LSF Stack plane"})
        kept.append(f)
    if not kept:
        raise RuntimeError("No frame could be registered to the existing LSF framing "
                           f"(min star pairs = {args.align_minpairs}).")
    msg = f"Alignment: {len(kept)}/{len(frames)} frame(s) registered to the LSF framing."
    if dropped:
        names = ", ".join(frames[i - 1]["filename"] or f"#{i}" for i in dropped)
        msg += f" Dropped (could not be registered, NOT stacked): {names}."
    log(msg, LogColor.SALMON if dropped else LogColor.GREEN)
    return kept


def _pick_ref_frame(frames: list, ref_frame: "int | None") -> int:
    """Batch-local index of the normalization reference: the user's --ref-frame (1-based
    sequence number) if given, else the best-quality frame."""
    if ref_frame is None:
        return _pick_best_quality_frame(frames)
    for n, f in enumerate(frames):
        if f.get("orig_index", f["index"]) == ref_frame - 1:
            return n
    raise RuntimeError(f"--ref-frame {ref_frame} is not among the included frames of the sequence.")


def _build_norm_baseline(args, frames: list, ref_idx: int, ref_pixels: np.ndarray,
                         channels: int, seqname: str) -> dict:
    """Per-channel (location, scale) of the baseline frame, stored in the file and used by every
    later addscale normalization. Written on initialization only."""
    loc, spr = [], []
    f = frames[ref_idx]
    for c in range(channels):
        plane = ref_pixels[..., c] if channels > 1 else ref_pixels
        st = f["stats_per_channel"][c] if c < len(f["stats_per_channel"]) else None
        est = lc.frame_estimators(plane, lc.valid_pixel_mask(plane), args.norm_estimator, st)
        if est is None or not est[1] > 0:
            raise RuntimeError(f"addscale: cannot estimate location/scale of the baseline frame "
                               f"{f['filename']} (channel {c}, estimator '{args.norm_estimator}').")
        loc.append(est[0])
        spr.append(est[1])
    return {"method": "addscale", "estimator": args.norm_estimator,
            "reference_frame": f"{seqname}:{f.get('orig_index', f['index'])}",
            "reference_filename": f.get("orig_filename", f["filename"]),
            "channels": channels, "loc": loc, "scale": spr}


def _load_norm_baseline(existing_data: dict, channels: int, args) -> dict:
    raw = existing_data["properties"].get("LSF:NormBaseline")
    if not raw:
        raise RuntimeError("Cannot append with norm_method='addscale': the existing file has no "
                           "LSF:NormBaseline block.")
    b = lc.loads(raw)
    if len(b.get("loc", [])) != channels or len(b.get("scale", [])) != channels:
        raise RuntimeError("Cannot append: the stored normalization baseline does not match the channel count.")
    if b.get("estimator") != args.norm_estimator:
        raise RuntimeError(f"Cannot append: baseline was measured with estimator '{b.get('estimator')}', "
                           f"this run uses '{args.norm_estimator}'.")
    return b


def _build_equalize_fits(args, frame: dict, ref_pixels: np.ndarray, channels: int) -> list:
    """Siril's -rgb_equal for the linearfit method, on initialization only.

    Siril (normalization.c, compute_factors_from_estimators) matches every channel's background
    location AND spread to the reference frame's GREEN channel:  s = scale_g / scale_c,
    pixel' = s * (pixel - loc_c) + loc_g.  That is what makes the per-channel noise of a Siril
    OSC stack come out about equal (the debayered green channel is intrinsically the quietest).

    Returns one NormFit (x*scale + offset) per channel: the map that equalizes the reference
    frame. On append nothing is needed: the accumulator the new frames are fitted against is
    already equalized."""
    est = []
    for c in range(channels):
        plane = ref_pixels[..., c]
        st = frame["stats_per_channel"][c] if c < len(frame["stats_per_channel"]) else None
        e = lc.frame_estimators(plane, lc.valid_pixel_mask(plane), args.norm_estimator, st)
        if e is None or not e[1] > 0:
            raise RuntimeError(f"equalize_rgb: cannot estimate location/scale of channel {c} of the "
                               f"reference frame {frame['filename']} (estimator '{args.norm_estimator}').")
        est.append(e)
    g_loc, g_sca = est[1]                                   # green, Siril's registration layer
    return [lc.addscale_fit(loc, sca, g_loc, g_sca) for loc, sca in est]


DEFAULT_SETTINGS = {
    "rejection": "sigma", "sigma_low": None, "sigma_high": None, "max_iters": 3,
    "gesdt_alpha": 0.05, "gesdt_max_frac": 0.3, "weighting": "none",
    "skip_normalization": False, "saturation_level": None, "output_norm": False,
    "align_transform": "homography", "align_interp": "lanczos4", "align_minpairs": 10,
    "norm_method": "linearfit", "norm_estimator": "lite",
    "equalize_rgb": None,      # None = auto: on for 3-channel (OSC) files, resolved in _run once channels are known
}


def _opt_float(text: str):
    text = (text or "").strip()
    return float(text) if text else None


def resolve_settings(args, stored: dict, log) -> None:
    """Fill every option the user did not give explicitly (None) from the
    existing file's stored settings, else from DEFAULT_SETTINGS. Explicit
    values always win, except weighting, which must match the file (mixed
    weighting schemes are not commensurate and would corrupt the accumulator)."""
    stored = stored or {}
    if (stored.get("weighting") and args.weighting is not None
            and args.weighting != stored["weighting"]):
        raise RuntimeError(
            f"Cannot append: file was built with weighting='{stored['weighting']}', "
            f"this run asks for '{args.weighting}'. Weights from different schemes are "
            f"not comparable.")
    if stored:      # appending: the normalization scheme is permanent (files without the key are linearfit)
        st_method = stored.get("norm_method") or "linearfit"
        if args.norm_method is not None and args.norm_method != st_method:
            raise RuntimeError(
                f"Cannot append: file was built with norm_method='{st_method}', this run asks for "
                f"'{args.norm_method}'. Coefficients from different methods are not commensurate.")
        if st_method == "addscale":
            if args.norm_estimator is not None and args.norm_estimator != stored.get("norm_estimator"):
                raise RuntimeError(f"Cannot append: norm_estimator is fixed at file creation "
                                   f"(file={stored.get('norm_estimator')}, now={args.norm_estimator}).")
        # RGB equalization changes the channels' relative scale and noise, so it is permanent for
        # every method (files written before it was recorded were not equalized).
        if args.equalize_rgb is not None and bool(args.equalize_rgb) != bool(stored.get("equalize_rgb", False)):
            raise RuntimeError(f"Cannot append: equalize_rgb is fixed at file creation "
                               f"(file={bool(stored.get('equalize_rgb', False))}, now={bool(args.equalize_rgb)}).")
    # sigma_low/high are only meaningful for the method they were stored with.
    same_method = args.rejection is None or args.rejection == stored.get("rejection")
    inherited, differs = [], []
    for key, default in DEFAULT_SETTINGS.items():
        cur = getattr(args, key, None)
        if cur is None:
            if stored.get(key) is not None and (same_method or key not in ("sigma_low", "sigma_high")):
                setattr(args, key, stored[key])
                inherited.append(f"{key}={stored[key]}")
            else:
                setattr(args, key, default)
        elif key in stored and stored[key] is not None and stored[key] != cur and key != "weighting":
            differs.append(f"{key}: file={stored[key]}, now={cur}")
    args.sigma_low, args.sigma_high = lc.resolve_rejection_params(
        args.rejection, args.sigma_low, args.sigma_high)
    if inherited:
        log("Settings inherited from existing LSF file: " + ", ".join(inherited))
    if differs:
        log("NOTE: settings differ from those the file was built with: " + "; ".join(differs),
            LogColor.SALMON)


class PhaseTimer:
    """Wall-clock accounting per pipeline phase, logged once at the end.

    Two ways to attribute time:
      add(name, dt)  time measured around an inner step (called many times in a loop);
      lap(name)      time since the previous lap, MINUS whatever add() already claimed
                     in that stretch, so sequential sections need no re-indentation and
                     nothing is counted twice.
    """

    def __init__(self):
        self.t0 = self._last = time.perf_counter()
        self._inner = 0.0
        self.totals = {}
        self.hidden = {}                 # worker-thread time that overlaps the main thread
        self._hlock = threading.Lock()

    def add(self, name: str, dt: float) -> None:
        self.totals[name] = self.totals.get(name, 0.0) + dt
        self._inner += dt

    def note(self, name: str, dt: float) -> None:
        """Time spent in a worker thread, overlapped with the main thread. Reported on its own
        (it is not additive with the table: it is the work the prefetch hides)."""
        with self._hlock:
            self.hidden[name] = self.hidden.get(name, 0.0) + dt

    def lap(self, name: str) -> None:
        now = time.perf_counter()
        dt = max(0.0, now - self._last - self._inner)
        self.totals[name] = self.totals.get(name, 0.0) + dt
        self._last, self._inner = now, 0.0

    def report(self, log) -> None:
        self.lap("finish")
        total = max(time.perf_counter() - self.t0, 1e-9)
        log(f"Timing: {total:.1f} s total")
        for name, dt in self.totals.items():
            if dt >= 0.05:
                log(f"  {name:<42s}{dt:9.1f} s  {100.0 * dt / total:5.1f}%")
        for name, dt in self.hidden.items():
            if dt >= 0.05:
                log(f"  (overlapped) {name:<29s}{dt:9.1f} s  in a worker thread, not in the total above")


class Prefetcher:
    """One-deep look-ahead. `load(i)` produces item i of n; while the caller works on item i,
    a single worker thread already loads item i+1, so reads overlap compute (NumPy and torch
    release the GIL in the heavy calls, and a blocked read is I/O wait anyway).

    get(i) must be called with i = 0, 1, 2, ... Each item is loaded exactly once, in order,
    by one thread at a time. At most two items are alive: the one being used and the next.
    An exception raised by `load` is re-raised by the get() of that item. With enabled=False
    nothing is threaded and get() just calls load() (identical results, no overlap).
    `on_time(seconds)` is called with the duration of every load (from the worker)."""

    def __init__(self, load, n: int, enabled: bool = True, on_time=None):
        self.load, self.n, self.on_time = load, n, on_time
        self.ex = ThreadPoolExecutor(max_workers=1, thread_name_prefix="lsf-prefetch") if enabled and n > 1 else None
        self.pending = {}
        self._sched(0)

    def _timed(self, i):
        t = time.perf_counter()
        try:
            return self.load(i)
        finally:
            if self.on_time is not None:
                self.on_time(time.perf_counter() - t)

    def _sched(self, i: int) -> None:
        if self.ex is not None and i < self.n and i not in self.pending:
            self.pending[i] = self.ex.submit(self._timed, i)

    def get(self, i: int):
        if self.ex is None:
            return self.load(i)
        self._sched(i)
        item = self.pending.pop(i).result()
        self._sched(i + 1)               # start the next load before the caller starts working
        return item

    def close(self) -> None:
        """Drop unstarted loads and wait for a running one (it may be using the Siril
        connection, which must be idle before the caller talks to Siril again)."""
        if self.ex is not None:
            self.ex.shutdown(wait=True, cancel_futures=True)
            self.pending.clear()


def warn_low_stack_count(siril: "s.SirilInterface", n_frames: int, threshold: int, log) -> bool:
    """FITS profile Section 10 (Initialize): if the first batch of a new file has fewer than
    `threshold` frames (default 25), warn the user before the file is created. The warning never
    blocks creation, and the profile defines no minimum for rejection, so rejection behaves the
    same at every count; the point is only that the per-pixel scatter it is based on is
    unreliable with so few frames, and a pixel rejected on that scatter stays rejected in the
    accumulator for good. Returns True when a warning was issued."""
    if not threshold or threshold <= 0 or n_frames >= threshold:
        return False
    msg = (f"Only {n_frames} frame(s) in the first batch of this new LSF file (fewer than {threshold}). "
           f"Per-pixel scatter estimated from so few frames is unreliable, and any rejection decision "
           f"made on it stays in the stack permanently. Consider starting the file with a larger "
           f"batch, or pass --min-frames-warn 0 to silence this warning. The file will still be created.")
    log("WARNING: " + msg, LogColor.SALMON)
    try:        # best effort: a visible, non-blocking pop-up when run from the Scripts menu
        if not siril.is_cli():
            siril.warning_messagebox(msg)
    except Exception:
        pass
    return True


def run(args: argparse.Namespace, siril: "s.SirilInterface") -> None:
    state = {"ctx": None, "timer": PhaseTimer()}
    ok = False
    try:
        _run(args, siril, state)
        ok = True
    finally:
        try:
            state["timer"].report(siril.log)
        except Exception:
            pass
        ctx = state.get("ctx")
        if ctx is not None:
            ctx.restore(siril, siril.log, ok)


def _run(args: argparse.Namespace, siril: "s.SirilInterface", state: dict) -> None:
    log = siril.log
    timer = state["timer"]

    backend = lc.configure_torch(getattr(args, "device", "auto"), getattr(args, "torch_precision", "auto"))
    lc.set_logger(log)
    log(f"Rejection backend: {backend}")

    # Resolve options first (alignment needs them): explicit > stored in --existing > defaults.
    if args.existing and args.append:
        if not os.path.exists(args.existing):
            raise RuntimeError(f"--existing file not found: {args.existing}")
        # Section 10, append step 1 / Section 8: validate before anything is merged. On a failed
        # check this raises with the check number, and the file is not touched.
        siril.update_progress("Validating existing LSF file", -1.0)
        lc.validate_lsf_file(args.existing)
        resolve_settings(args, lc.read_lsf_settings(args.existing)["settings"], log)
    else:
        resolve_settings(args, {}, log)

    if args.seq:
        siril.cmd("load_seq", args.seq)
    seq = siril.get_seq()
    if seq is None:
        raise RuntimeError("No sequence is loaded in Siril (pass --seq or load one first).")
    seqname = args.seq or seq.seqname

    included = [i for i in range(seq.number) if seq.imgparam[i].incl]
    if not included:
        raise RuntimeError(f"Sequence '{seqname}' has no included frames to stack.")

    channels = max(1, int(getattr(seq, "nb_layers", 1)))

    log(f"LSF stack: '{seqname}', {len(included)}/{seq.number} frames included, "
        f"{channels} channel(s).")
    siril.update_progress("Reading frames", -1.0)

    # Metadata only; pixels are streamed one frame at a time later.
    frames = []
    for n, idx in enumerate(included):
        frames.append(gather_frame_meta(siril, idx, channels))
        siril.update_progress(f"Reading frame metadata {n + 1}/{len(included)}",
                              (n + 1) / len(included) * 0.1)

    timer.lap("Siril sequence + frame metadata")

    if args.existing and args.append and not args.skip_align:
        frames = align_to_lsf(siril, args, frames, seqname, channels, state, log)
    elif args.existing and args.append:
        log("Alignment skipped (--skip-align): frames are added as they are.", LogColor.SALMON)

    timer.lap("alignment (Siril)")

    # Frame geometry from the sequence that is loaded NOW (after alignment this is the
    # registered temp sequence, which has the LSF grid's size). Falls back to reading
    # the first frame only when the sequence does not report usable dimensions or has
    # variable frame sizes.
    cur_seq = siril.get_seq() or seq
    rx, ry = int(getattr(cur_seq, "rx", 0) or 0), int(getattr(cur_seq, "ry", 0) or 0)
    if rx > 0 and ry > 0 and not getattr(cur_seq, "is_variable", False):
        width, height = rx, ry
        actual_channels = max(1, int(getattr(cur_seq, "nb_layers", channels) or channels))
    else:
        probe = read_frame_pixels(siril, frames[0]["index"])
        height, width = probe.shape[0], probe.shape[1]
        actual_channels = probe.shape[2] if probe.ndim == 3 else 1
        del probe
    if actual_channels != channels:
        log(f"WARNING: sequence reports nb_layers={channels} but the first frame's pixel "
            f"data has {actual_channels} channel(s); trusting the pixel data.", LogColor.SALMON)
        channels = actual_channels
        for f in frames:   # per-channel stats were gathered with the old count
            f["stats_per_channel"] = [None] * channels if channels != len(f["stats_per_channel"]) else f["stats_per_channel"]

    if args.equalize_rgb is None:       # auto: like Siril's OSC scripts (-rgb_equal); never on an old file
        args.equalize_rgb = bool(channels == 3 and not (args.existing and args.append))
    elif args.equalize_rgb and channels != 3:
        log("equalize_rgb ignored: the sequence is not 3-channel.", LogColor.SALMON)
        args.equalize_rgb = False
    if channels == 3:
        log("RGB equalization to green (Siril -rgb_equal): " + ("on" if args.equalize_rgb else "off")
            + ("" if args.equalize_rgb or not (args.existing and args.append)
               else " (as the file was created)") + ".")

    sensor_profile = _infer_sensor_profile(channels, frames[0]["header"], args)
    filter_name = _extract_filter(frames[0]["header"])
    canonical_grid = {
        "width": width,
        "height": height,
        "channels": channels,
        "wcs": _extract_wcs(frames[0]["header"]),
        "filter": filter_name,
        **sensor_profile,
    }
    log(f"Sensor profile: {sensor_profile['sensor_class']}"
        + (f" ({sensor_profile['bayer_pattern']})" if sensor_profile["sensor_class"] == "cfa" else "")
        + f", demosaic={sensor_profile['demosaic_method']}, filter={filter_name or '(none)'}.")

    # ---- Load / validate existing accumulator, or start fresh ------------
    existing_data = None
    prev_meta = None
    if args.existing and args.append:
        existing_data = lc.read_lsf(args.existing)
        prev_meta = existing_data.get("meta") or None      # LSFUID, LSFGEN, unknown content (C7)
        lc.check_lsf_file(existing_data, args.existing)
        ex_grid = lc.loads(existing_data["properties"]["LSF:CanonicalGrid"])
        if ex_grid["width"] != width or ex_grid["height"] != height:
            raise RuntimeError(
                f"Cannot append: existing LSF grid is {ex_grid['width']}x{ex_grid['height']}, "
                f"this batch is {width}x{height}."
            )
        if ex_grid["channels"] != channels:
            raise RuntimeError(
                f"Cannot append: existing LSF file has {ex_grid['channels']} channel(s), "
                f"this batch has {channels}."
            )
        # Section 9 / 10.1: sensor class, CFA pattern, demosaic method and
        # filter are permanent once the file is created -- hard errors,
        # not warnings, on mismatch.
        ex_sensor = ex_grid["sensor_class"]
        if ex_sensor != sensor_profile["sensor_class"]:
            raise RuntimeError(
                f"Cannot append: existing LSF file is sensor_class='{ex_sensor}', "
                f"this batch is '{sensor_profile['sensor_class']}'."
            )
        if ex_sensor == "cfa":      # (this batch is cfa too: a mismatch raised above)
            ex_bayer = ex_grid["bayer_pattern"]
            new_bayer = sensor_profile["bayer_pattern"]
            if "unknown" not in (ex_bayer, new_bayer) and ex_bayer != new_bayer:
                raise RuntimeError(
                    f"Cannot append: existing LSF file is Bayer pattern '{ex_bayer}', "
                    f"this batch is '{new_bayer}'."
                )
        ex_demosaic = ex_grid["demosaic_method"]
        if ex_demosaic != sensor_profile["demosaic_method"]:
            raise RuntimeError(
                f"Cannot append: existing LSF file used demosaic_method='{ex_demosaic}' "
                f"(Section 10.3 -- fixed at file creation), this batch is "
                f"'{sensor_profile['demosaic_method']}'."
            )
        ex_filter = ex_grid["filter"]      # None when the frames carried no FILTER keyword
        if ex_filter and filter_name and ex_filter != filter_name:
            raise RuntimeError(
                f"Cannot append: existing LSF file is filter='{ex_filter}' "
                f"(Section 10.1 -- a file represents one filter), this batch is "
                f"filter='{filter_name}'."
            )
        old_scale = _wcs_pixel_scale(ex_grid["wcs"])
        new_scale = _wcs_pixel_scale(canonical_grid["wcs"])
        if old_scale and new_scale and old_scale > 0:
            rel_diff = abs(old_scale - new_scale) / old_scale
            if rel_diff > 0.01:
                log(f"WARNING: pixel scale differs from existing LSF file by {rel_diff*100:.2f}% "
                    f"-- grids may not actually be compatible.", LogColor.SALMON)
        log(f"Appending to existing LSF file: {args.existing}")
    elif args.existing and not args.append:
        log("--existing given without --append: ignoring its history, "
            "this will be an independent, fresh stack.", LogColor.SALMON)

    # ---- Initialization only: low stack count warning (profile Section 10) ----
    if existing_data is None:
        warn_low_stack_count(siril, len(frames), getattr(args, "min_frames_warn", LOW_STACK_WARN_FRAMES), log)

    # ---- Output-norm transform and canonical grid (append) ----
    was_append = existing_data is not None
    old_norm = None            # (mini, maxi) if the existing Stack plane is normalized
    if existing_data is not None:
        info = lc.loads(existing_data["properties"]["LSF:OutputNorm"])
        if info.get("enabled"):
            old_norm = (float(info["mini"]), float(info["maxi"]))
        flagged = any(k == "LSFNORM" and v is True
                      for k, v, _ in existing_data["images"]["Stack"]["fits_keywords"])
        if flagged and old_norm is None:
            raise RuntimeError(
                "Cannot append: the file's Stack plane is flagged as normalized (LSFNORM) but the "
                "LSF:OutputNorm transform is missing, so the raw accumulator mean cannot be "
                "recovered. Refusing to merge into a mis-scaled mean.")
        # Spec Section 7: the canonical grid is fixed at creation. Carry it forward.
        canonical_grid = ex_grid
    if old_norm is not None:
        log(f"Existing Stack plane is output-normalized (min={old_norm[0]:.6g}, max={old_norm[1]:.6g}); "
            f"inverting for merge.")

    # ---- Ingestion-time normalization reference (Section 8.3) -------------
    # Append: the accumulator's own current best estimate (sum/weight),
    # restricted to pixels it has actually seen (weight > 0). Fresh init:
    # the best-quality frame in this batch, chosen via Section 8.1 metrics
    # -- that frame's own scale becomes this file's permanent baseline.
    # NOTE: These reference arrays are used only during Pass 1 (ingestion-time
    # normalization) and freed before Pass 2 to reduce peak memory usage.
    timer.lap("open/validate existing file")
    old_ledger = lc.loads(existing_data["properties"]["LSF:FrameLedger"]) if existing_data is not None else []
    reference_planes = []       # per channel: 2D array
    reference_valid_masks = []  # per channel: 2D bool array
    reference_is_accumulator = bool(old_ledger)       # every ledgered frame contributed to the accumulator

    norm_baseline = None
    eq_fits = None              # linearfit + equalize_rgb, initialization only: per-channel equalizing map
    if reference_is_accumulator and args.norm_method == "addscale":
        norm_baseline = _load_norm_baseline(existing_data, channels, args)
        if getattr(args, "ref_frame", None):
            log("--ref-frame ignored: the stored addscale baseline is used.", LogColor.SALMON)
        if args.norm_estimator == "siril" and not args.skip_align:
            log("NOTE: estimator 'siril' uses statistics measured BEFORE alignment resampling; "
                "use --skip-align or --norm-estimator lite for exact parity.", LogColor.SALMON)
        log(f"Normalization reference: stored addscale baseline (frame {norm_baseline.get('reference_frame')}).")
    elif reference_is_accumulator:
        stack_data = existing_data["images"]["Stack"]["data"]
        weight_data = existing_data["images"]["WeightSum"]["data"]
        for c in range(channels):
            plane = stack_data[..., c] if channels > 1 else stack_data
            wplane = weight_data[..., c] if channels > 1 else weight_data
            reference_planes.append(plane)
            reference_valid_masks.append(wplane > 0)
        log("Normalization reference: existing accumulator (sum/weight).")
    else:
        ref_idx = _pick_ref_frame(frames, getattr(args, "ref_frame", None))
        ref_pixels = read_frame_pixels(siril, frames[ref_idx]["index"])
        for c in range(channels):
            plane = ref_pixels[..., c] if channels > 1 else ref_pixels
            reference_planes.append(plane)
            reference_valid_masks.append(lc.valid_pixel_mask(plane, args.saturation_level))
        if args.norm_method == "addscale":
            norm_baseline = _build_norm_baseline(args, frames, ref_idx, ref_pixels, channels, seqname)
            log("addscale baseline (per channel) loc=" + ", ".join(f"{v:.6g}" for v in norm_baseline["loc"])
                + " scale=" + ", ".join(f"{v:.6g}" for v in norm_baseline["scale"])
                + f" [{args.norm_estimator}]" + (", equalized to green" if args.equalize_rgb and channels == 3 else ""))
        if (args.equalize_rgb and channels == 3 and args.norm_method == "linearfit"
                and not args.skip_normalization):
            eq_fits = _build_equalize_fits(args, frames[ref_idx], ref_pixels, channels)
            log("RGB equalization (linearfit): reference channel scales vs green = "
                + ", ".join(f"{f.scale:.4g}" for f in eq_fits))
        log(f"Normalization reference: batch frame #{ref_idx} ({frames[ref_idx]['filename']}), "
            + ("chosen with --ref-frame." if getattr(args, "ref_frame", None)
               else "chosen as the best-quality frame in this batch."))

    timer.lap("normalization reference setup")

    # ---- Pass 1: normalize one frame at a time into a disk-backed spill ----
    # The spill lives in its own throw-away folder (--temp-dir, default
    # <output dir>/lsf_temp) and is ALWAYS removed in the finally block below.
    siril.update_progress("Normalizing frames", 0.15)
    n_frames = len(frames)
    norm_shape = (n_frames, height, width, channels) if channels > 1 else (n_frames, height, width)
    plane_shape = (height, width, channels) if channels > 1 else (height, width)
    frame_weights = np.empty((n_frames, channels), dtype=np.float64)
    ledger_rows = []

    def norm_ref_channel(c: int) -> int:
        # Siril -rgb_equal: every channel is matched to the registration layer (green, 1).
        return 1 if (args.equalize_rgb and channels == 3) else c

    default_tmp = os.path.join(os.path.dirname(os.path.abspath(args.output)), "lsf_temp")
    spill_root = args.temp_dir or default_tmp
    created_spill_root = not os.path.isdir(spill_root)
    os.makedirs(spill_root, exist_ok=True)
    spill_dir = tempfile.mkdtemp(prefix="lsf_spill_", dir=spill_root)
    spill_path = os.path.join(spill_dir, "normalized.f32")
    log(f"Temporary normalized-frame spill: {spill_path} "
        f"(~{np.prod(norm_shape) * 4 / 2**30:.1f} GB, removed when done)")

    normalized = None
    pf1 = pf2 = None
    prefetch = getattr(args, "prefetch", "spill")
    try:
        normalized = np.memmap(spill_path, dtype=np.float32, mode="w+", shape=norm_shape)
        pf1 = Prefetcher(
            lambda i: read_frame_pixels(siril, frames[i]["index"], frames[i].get("pixel_scale", 1.0)),
            n_frames, enabled=(prefetch == "all"),
            on_time=lambda dt: timer.note("pass 1: Siril pixel read (worker)", dt))

        for n, f in enumerate(frames):
            _t = time.perf_counter()
            px = pf1.get(n)          # with --prefetch all this only waits for the look-ahead read
            timer.add("pass 1: Siril pixel read", time.perf_counter() - _t)
            px_channels = px.shape[2] if px.ndim == 3 else 1
            if px_channels != channels:
                raise RuntimeError(
                    f"Frame {f['filename']} has {px_channels} channel(s), expected {channels}.")
            if px.shape[0] != height or px.shape[1] != width:
                raise RuntimeError(
                    f"Frame {f['filename']} has shape {px.shape[:2]}, expected "
                    f"({height}, {width}). All frames must already share the canonical grid "
                    f"(register/resample before calling this script)."
                )
            norm_scales, norm_offsets, norm_methods = [], [], []
            norm_loc, norm_spread = [], []
            for c in range(channels):
                src_c = px[..., c] if channels > 1 else px
                # Covered pixels only (finite and non-zero). Empty border
                # pixels become NaN so the batch stats ignore them.
                src_mask_c = lc.valid_pixel_mask(src_c)
                norm_loc.append(None)
                norm_spread.append(None)
                if args.skip_normalization:
                    scale, offset, method = 1.0, 0.0, "skipped"
                elif args.norm_method == "addscale":
                    _t = time.perf_counter()
                    rc = norm_ref_channel(c)
                    st_c = f["stats_per_channel"][c] if c < len(f["stats_per_channel"]) else None
                    est = lc.frame_estimators(src_c, src_mask_c, args.norm_estimator, st_c)
                    if est is None:
                        scale, offset, method = 1.0, 0.0, "failed"
                        log(f"  WARNING: {f['filename']} channel {c}: no '{args.norm_estimator}' "
                            f"location/scale estimate -- left unnormalized.", LogColor.SALMON)
                    else:
                        fit = lc.addscale_fit(est[0], est[1], norm_baseline["loc"][rc], norm_baseline["scale"][rc])
                        scale, offset, method = fit.scale, fit.offset, fit.method
                        norm_loc[-1], norm_spread[-1] = est
                    timer.add("pass 1: normalization fit", time.perf_counter() - _t)
                else:
                    src_mask_fit = (src_mask_c if args.saturation_level is None
                                    else src_mask_c & (src_c < args.saturation_level))
                    _t = time.perf_counter()
                    # Reference noise: the accumulator mean is far cleaner than one frame
                    # (plain least squares is enough); on first init the reference is a single
                    # batch frame, so both axes are noisy and the fit is errors-in-variables.
                    # If the stored Stack plane is output-normalized, the fit converts the
                    # reference to raw units itself, so scale/offset come out in raw units.
                    ref_to_raw = ((old_norm[1] - old_norm[0], old_norm[0])
                                  if (reference_is_accumulator and old_norm is not None) else None)
                    fit = lc.fit_normalization(
                        src_c, reference_planes[c], reference_valid_masks[c] & src_mask_fit,
                        reference_is_clean=reference_is_accumulator, reference_to_raw=ref_to_raw)
                    timer.add("pass 1: normalization fit", time.perf_counter() - _t)
                    scale, offset, method = fit.scale, fit.offset, fit.method
                    if eq_fits is not None:
                        # frame -> reference (scale, offset), then reference -> equalized reference
                        e = eq_fits[c]
                        scale, offset = e.scale * scale, e.scale * offset + e.offset
                    if method == "offset_only":
                        log(f"  NOTE: {f['filename']} channel {c}: {fit.reason}; matched the sky "
                            f"level only (scale=1, offset={offset:.6g}). Any transparency change "
                            f"between this frame and the reference is not corrected.", LogColor.SALMON)
                    elif method == "failed":
                        log(f"  WARNING: {f['filename']} channel {c}: {fit.reason} -- left "
                            f"unnormalized (scale=1, offset=0) for this channel.", LogColor.SALMON)
                norm_methods.append(method)
                norm_scales.append(scale)
                norm_offsets.append(offset)
                _t = time.perf_counter()
                normalized_c = np.where(
                    src_mask_c,
                    lc.apply_linear_normalization(src_c, scale, offset),
                    np.nan).astype(np.float32)
                if channels > 1:
                    normalized[n, ..., c] = normalized_c
                else:
                    normalized[n] = normalized_c
                timer.add("pass 1: apply + spill write", time.perf_counter() - _t)
            del px
            siril.update_progress(f"Normalizing frame {n + 1}/{n_frames}",
                                  0.15 + 0.45 * (n + 1) / n_frames)

            reg = f["reg"]
            stats_pc = f["stats_per_channel"]
            nbstars = getattr(reg, "number_of_stars", None) if reg else None
            wfwhm = getattr(reg, "weighted_fwhm", None) if reg else None
            weights_this_frame, background_medians = [], []
            noises, noises_norm = [], []
            frame_stats = {k: [] for k in LEDGER_STAT_FIELDS}
            for c in range(channels):
                st = stats_pc[c] if c < len(stats_pc) else None
                noise_c = getattr(st, "bgnoise", None) if st else None
                for k in LEDGER_STAT_FIELDS:
                    frame_stats[k].append(_num(getattr(st, k, None)) if st else None)
                _n = _num(noise_c)
                noises.append(_n)
                noises_norm.append(None if _n is None else _n * abs(norm_scales[c]))
                # bgnoise was measured before normalization; the accumulator sees the data
                # after `pixel*scale + offset`, so weight by the noise in those units.
                weights_this_frame.append(lc.compute_frame_weight(
                    args.weighting, noise_c, nbstars, wfwhm, scale=norm_scales[c]))
                background_medians.append(getattr(st, "median", None) if st else None)
            frame_weights[n] = weights_this_frame

            header = f["header"]
            homography = f.get("alignment_H")        # transform onto the LSF grid wins
            if homography is None:
                homography = _h_list(reg)

            ledger_rows.append({
                "frame_id": f"{seqname}:{f.get('orig_index', f['index'])}",
                "original_filename": f.get("orig_filename", f["filename"]),
                "capture_timestamp": header.get("DATE-OBS"),
                "exposure_seconds": header.get("EXPTIME", header.get("EXPOSURE")),
                "gain": header.get("GAIN"),
                "offset": header.get("OFFSET"),
                "sensor_temp_c": header.get("CCD-TEMP"),
                "filter": header.get("FILTER"),
                "original_header": f["header_str"],
                "registration_transform": homography,
                "alignment": f.get("alignment"),
                "normalization": dict(
                    {"scale": norm_scales, "offset": norm_offsets, "method": norm_methods},
                    **({"estimator": {"kind": args.norm_estimator, "loc": norm_loc, "spread": norm_spread}}
                       if args.norm_method == "addscale" else {})),
                "weight_applied": weights_this_frame,
                "quality": {
                    "fwhm": _num(getattr(reg, "fwhm", None)) if reg else None,
                    "weighted_fwhm": _num(wfwhm),
                    "roundness": _num(getattr(reg, "roundness", None)) if reg else None,
                    "star_count": getattr(reg, "number_of_stars", None) if reg else None,
                    "registration_quality": _num(getattr(reg, "quality", None)) if reg else None,
                    "registration_background": _num(getattr(reg, "background_lvl", None)) if reg else None,
                    # Per-plane lists (length 1 mono, 3 OSC), as for the other plane fields.
                    "background_median": background_medians,
                    # Siril background noise (bgnoise) as measured on the frame BEFORE
                    # normalization, and the same value expressed in the accumulator's
                    # units (x |normalization scale|). The latter is what an inverse-variance
                    # weight would need: w = 1 / noise_normalized**2.
                    "noise": noises,
                    "noise_normalized": noises_norm,
                    "stats": frame_stats,
                },
            })

        _t = time.perf_counter()
        normalized.flush()
        timer.add("pass 1: spill flush", time.perf_counter() - _t)
        pf1.close()
        timer.lap("pass 1 other (stats, ledger rows)")

        # Free normalization reference arrays before Pass 2 to reduce peak memory.
        # These were only needed for ingestion-time normalization (Pass 1, line 639).
        # They are NOT needed for Siril-style output normalization, which operates
        # only on the final stacked result.
        reference_planes = []
        reference_valid_masks = []
        gc.collect()

        # ---- Pass 2: rejection + combine in row strips ------------------------
        # Rejection/weighting are per-pixel across frames, so each row band is
        # independent: only (N, rows, W[, C]) is ever resident, never the full stack.
        # float64 temporaries inside compute_batch_stats; the torch path keeps more
        # (sort indices, masks) alive at once, so budget conservatively.
        BYTES_PER_SAMPLE = lc.bytes_per_sample(args.rejection)
        if prefetch != "off":
            # one extra strip alive (4 B/sample) plus the prefetched slice of the old planes
            # (18 B per pixel-channel: mean, M2, W, W2 float32 + count uint16, spread over the N
            # frames of each pixel)
            BYTES_PER_SAMPLE += 4 + ((18 // max(1, n_frames) + 1) if existing_data is not None else 0)
        budget = int(args.mem_budget_gb * 2 ** 30)
        rows = max(1, min(height, budget // max(1, n_frames * width * channels * BYTES_PER_SAMPLE)))
        log(f"Rejecting & combining in strips of {rows} row(s) "
            f"({n_frames} frames, ~{args.mem_budget_gb:.1f} GB working budget).")

        final_mean = np.empty(plane_shape, dtype=np.float32)
        final_m2 = np.empty(plane_shape, dtype=np.float32)
        final_w = np.empty(plane_shape, dtype=np.float32)
        final_w2 = np.empty(plane_shape, dtype=np.float32)
        final_count = np.empty(plane_shape, dtype=np.uint16)

        if existing_data is not None:
            old_mean = existing_data["images"]["Stack"]["data"]
            old_weight = existing_data["images"]["WeightSum"]["data"]
            old_m2 = existing_data["images"]["SumSq"]["data"]
            old_w2 = existing_data["images"]["SumW2"]["data"]
            old_count = existing_data["images"]["ContributionCount"]["data"]

        starts = list(range(0, height, rows))

        def load_band(i):
            """Everything strip i needs from disk, copied into RAM: the normalized frames and,
            on append, the matching rows of the four old accumulator planes."""
            a = starts[i]
            b = min(height, a + rows)
            band_ = np.array(normalized[:, a:b])    # contiguous read, copied into RAM
            old_ = None
            if existing_data is not None:
                # lc.native: a FITS plane is a big-endian memory map; the strip is copied into RAM
                # in native byte order (torch and the merge arithmetic want that).
                old_ = tuple(lc.native(p[a:b]) for p in (old_mean, old_m2, old_weight, old_w2, old_count))
            return band_, old_

        pf2 = Prefetcher(load_band, len(starts), enabled=(prefetch != "off"),
                         on_time=lambda dt: timer.note("pass 2: spill + old-plane read (worker)", dt))
        for bi, r0 in enumerate(starts):
            r1 = min(height, r0 + rows)
            _t = time.perf_counter()
            band, old_b = pf2.get(bi)   # with prefetch this only waits for the look-ahead read
            timer.add("pass 2: spill + old-plane read", time.perf_counter() - _t)
            for c in range(channels):
                strip = band[..., c] if channels > 1 else band
                _t = time.perf_counter()
                bm, bm2, bw, bw2, bc = lc.compute_batch_stats(
                    strip, frame_weights[:, c], rejection=args.rejection,
                    sigma_low=args.sigma_low, sigma_high=args.sigma_high,
                    max_iters=args.max_iters, gesdt_alpha=args.gesdt_alpha,
                    gesdt_max_frac=args.gesdt_max_frac)
                timer.add("pass 2: rejection + batch stats", time.perf_counter() - _t)
                _t = time.perf_counter()
                osl = (slice(r0, r1), slice(None)) + ((c,) if channels > 1 else ())
                if existing_data is not None:
                    lsl = (slice(None), slice(None)) + ((c,) if channels > 1 else ())   # within the strip
                    om = old_b[0][lsl]
                    if old_norm is not None:
                        om = lc.invert_output_norm(om, old_b[2][lsl], old_norm[0], old_norm[1])
                    bm, bm2, bw, bw2, bc = lc.combine_weighted_stats(
                        om, old_b[1][lsl], old_b[2][lsl], old_b[3][lsl], old_b[4][lsl],
                        bm, bm2, bw, bw2, bc)
                timer.add("pass 2: merge with history", time.perf_counter() - _t)
                final_mean[osl], final_m2[osl] = bm, bm2
                final_w[osl], final_w2[osl], final_count[osl] = bw, bw2, bc
            del band, old_b
            siril.update_progress("Rejecting & combining", 0.6 + 0.3 * r1 / height)
        timer.lap("pass 2 other (strip bookkeeping)")
        _ts = lc.torch_summary()
        if _ts:
            log(_ts, LogColor.SALMON)
    finally:
        # Stop the look-ahead workers first (one may be mid-read), then remove the temporary
        # spill, on success AND on error.
        for _pf in (pf1, pf2):
            if _pf is not None:
                _pf.close()
        normalized = None
        gc.collect()
        shutil.rmtree(spill_dir, ignore_errors=True)
        if created_spill_root:
            try:
                os.rmdir(spill_root)   # only succeeds if empty
            except OSError:
                pass

    timer.lap("spill cleanup")
    full_ledger = old_ledger + ledger_rows

    norm_info = {"enabled": False}
    if args.output_norm:
        # Siril-style min-max normalization of the *stored* Stack plane only.
        # SumSq/WeightSum stay raw; (mini, maxi) is recorded so --append can invert it.
        res = lc.apply_output_norm(final_mean, final_w)
        if res is None:
            log("Output normalization skipped: no valid dynamic range.", LogColor.SALMON)
        else:
            norm_info = {"enabled": True, "mini": res[0], "maxi": res[1]}

    timer.lap("output normalization")

    # Release every memory map onto the existing file before (possibly) overwriting it
    # in place (--output == --existing); required on Windows. Setting existing_data = None
    # alone is NOT enough: the old_* arrays, the reference-block locals and the last strip's
    # `om` view all hold views of the same maps (the reference lists were already released
    # before Pass 2), so every one of them must be dropped before the mapping is closed.
    old_mean = old_weight = old_m2 = old_w2 = old_count = None
    stack_data = weight_data = plane = wplane = om = None
    existing_data = None
    gc.collect()

    siril.update_progress("Writing LSF file", 0.9)

    # ---- Assemble and write the output file --------------------------------
    total_exposure = 0.0
    for row in full_ledger:
        try:
            total_exposure += float(row.get("exposure_seconds") or 0.0)
        except (TypeError, ValueError):
            pass
    dates = [r["capture_timestamp"] for r in full_ledger if r.get("capture_timestamp")]

    fits_kw = []
    wcs = canonical_grid["wcs"]
    for k, v in wcs.items():
        fits_kw.append((k, v, "canonical grid WCS"))
    fits_kw.append(("STACKCNT", len(full_ledger), "total frames ever contributed"))
    fits_kw.append(("EXPTIME", total_exposure, "total exposure, seconds"))
    if dates:
        fits_kw.append(("DATE-OBS", min(dates), "earliest contributing frame"))
    if args.object_name:
        fits_kw.append(("OBJECT", args.object_name, ""))
    if filter_name:
        fits_kw.append(("FILTER", filter_name, ""))
    # No LSF* keywords and no BAYERPAT here: the container derives LSF / LSFVER / LSFUID / LSFGEN / LSFSENS /
    # LSFBAYR / LSFDMOD / LSFNORM / LSFNMIN / LSFNMAX / LSFDSUM from the properties below (FITS profile,
    # Section 4: an OSC primary is already RGB, so the pattern lives in LSFBAYR only).
    fits_kw.append(("HISTORY", f"lsf_stack.py: {'appended to' if was_append else 'initialized from'} "
                                f"'{seqname}' ({len(frames)} frames) at {time.strftime('%Y-%m-%dT%H:%M:%S')}", ""))

    images = [
        lc.LsfImage("Stack", final_mean, fits_keywords=fits_kw),
        lc.LsfImage("WeightSum", final_w),
        lc.LsfImage("SumSq", final_m2),
        lc.LsfImage("SumW2", final_w2),
        lc.LsfImage("ContributionCount", final_count),
    ]
    properties = {
        "LSF:FrameLedger": lc.dumps(full_ledger),
        "LSF:CanonicalGrid": lc.dumps(canonical_grid),
        "LSF:Settings": lc.dumps({k: getattr(args, k) for k in DEFAULT_SETTINGS}),
        "LSF:OutputNorm": lc.dumps(norm_info),
    }
    if norm_baseline is not None:
        properties["LSF:NormBaseline"] = lc.dumps(norm_baseline)

    timer.lap("assemble header + ledger JSON")
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    # Write to a temp file in the same folder, then atomically swap it in, so a failure
    # mid-write can never truncate the only copy when --output == --existing.
    tmp_out = args.output + ".lsf_tmp"
    # ROWORDER: a new file takes the source frames' (else the profile default); an append keeps the file's.
    row_order = None if was_append else (frames[0]["header"].get("ROWORDER") or None)
    try:
        wstats = lc.write_lsf(tmp_out, images, properties, compress="gzip2" if args.compress else None,
                              meta=prev_meta, new_rows=len(ledger_rows), row_order=row_order)
        if args.keep_backup and os.path.exists(args.output):      # Section 10: offer to keep the previous file
            os.replace(args.output, args.output + ".bak")
        os.replace(tmp_out, args.output)                          # atomic swap (the temp file was fsync'ed)
        timer.lap("write LSF file")
    finally:
        if os.path.exists(tmp_out):
            try:
                os.remove(tmp_out)
            except OSError:
                pass

    siril.update_progress("Done", 1.0)
    log(f"LSF stack written to {args.output} "
        f"({len(full_ledger)} total frames in history, {len(frames)} added this run; "
        f"{wstats['bytes'] / 2**20:.1f} MB, {'compressed ' + wstats['mode'] if wstats['mode'] else 'uncompressed'}"
        + (f", LSFGEN {wstats['gen']}" if wstats.get("gen") else "") + ").",
        LogColor.GREEN)


def main():
    siril = s.SirilInterface()
    siril.connect()
    try:
        parser = build_arg_parser()
        argv = sys.argv[1:]
        current_seqname = None
        try:
            seq = siril.get_seq()
            if seq:
                current_seqname = seq.seqname
        except Exception:
            pass

        if not argv and not siril.is_cli():
            args = run_gui_and_get_args(current_seqname)
            if args is None:
                siril.log("LSF stack cancelled by user.")
                return
        else:
            args = parser.parse_args(argv)
            if not args.seq:
                args.seq = current_seqname
            if not args.output:
                raise RuntimeError("--output is required (or use the GUI form).")

        run(args, siril)
    except Exception as e:
        siril.log(f"LSF stack failed: {e}", LogColor.RED)
        traceback.print_exc()
        try:
            siril.error_messagebox(f"LSF stack failed:\n{e}")
        except Exception:
            pass
        sys.exit(1)
    finally:
        siril.disconnect()


if __name__ == "__main__":
    main()