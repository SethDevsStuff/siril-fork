#!/usr/bin/env python3
"""
lsf_compare.py -- end-to-end test of Siril's native Living Stack Format (LSF v0.4) commands
(lsf_init, lsf_append, lsf_export): fresh LSF stacks versus iterative appends, parity against
Siril's own `stack` command, v0.4 structural checks, SSIM / MS-SSIM between the full stacks, and a
noise estimate for every product.

Author: (fill in your name here before distributing)
License: (fill in a license here before distributing, e.g. GPL-3.0)
Contact: (fill in a contact / forum / repo URL here before distributing)

INPUT
-----
A Siril sequence whose frames are already calibrated and registered to a common framing (this
script does no registration; the appends use lsf_append -noreg). Only the frames SELECTED
(included) in the sequence are used. They are split into 4 sets of nearly equal size
(sizes differ by at most one frame; contiguous blocks by default, or --split interleaved).
Every stacking command works on the current selection: lsf_init / lsf_append use the included
frames, `stack` gets -filter-included.
CFA drizzle files are not covered: lsf_append -noreg refuses them, and lsf_init refuses a CFA or
drizzled sequence without -drizzle, which this script does not pass.

PRODUCTS  (all written to --outdir, default <Siril working dir>/lsf_compare; prefix --prefix)
--------
Stage base   (lsf_init, fresh stacks)
             <p>_lsf_set_1.lsf.fits ... <p>_lsf_set_4.lsf.fits   one per set
             <p>_lsf_cum_1/12/123/1234.lsf.fits                   running stacks over sets 1..k
Stage siril  (native `stack ... -32b -filter-included`, same settings)
             <p>_siril_set_1.fit ... <p>_siril_set_4.fit
             <p>_siril_cum_1/12/123/1234.fit
Stage sos    <p>_lsf_sos.lsf.fits   lsf_init on a scratch sequence made of the four per-set
                                    stacks (exported with lsf_export)
Stage iter   <p>_lsf_1.lsf.fits     = set 1 (copy of lsf_set_1)
             <p>_lsf_2..4.lsf.fits  previous file copied, then `lsf_append <copy> <seq> -noreg`
                                    with the next set selected (lsf_append writes in place)
Stage check  <p>_check.csv          pure-Python v0.4 structure checks of every .lsf.fits product:
                                    EXTNAMEs (LSF_ACC, LSF_COUNT, LSF_SETTINGS, LSF_LEDGER last),
                                    LSF_ACC NAXIS3 == 4*P (5*P with the EXP planes of an exposure
                                    map), LSFUID/LSFGEN identical across LSF HDUs,
                                    LSFGEN +1 along lsf_1..lsf_4, STACKCNT == ledger rows,
                                    LSFHASH == sha256(primary data), plus `fitsverify` if on PATH.
                                    Raises if any check fails.
Stage parity <p>_parity.csv         lsf_set_k vs siril_set_k and lsf_cum_X vs siril_cum_X:
                                    per-channel max abs and max rel difference; pixels that
                                    differ at full and at partial coverage (LSF_COUNT vs
                                    STACKCNT) and that have data in only one image; Siril's
                                    Noise Estimation of both images and the rows where it
                                    differs. These are the same frames with the same settings,
                                    so the images are expected to be identical. Known, expected
                                    exceptions (logged): wfwhm / nbstars weighting, and sigmedian
                                    rejection without weights.
             <p>_diff_<pair>.fit    LSF - Siril, written only for pairs that differ.
Stage ssim   <p>_ssim.csv           SSIM and MS-SSIM, per channel, between every pair of the four full
                                    stacks: lsf_cum_1234 (fresh), lsf_4 (iterative), lsf_sos (stack on
                                    stack) and siril_cum_1234. See SSIM below.
             <p>_ssim_<pair>.fit    1 - SSIM map, only with --ssim-maps.
Stage flip  <p>_flip.csv           On a scratch copy of an LSF product: its
                                    Stack flipped with the old DATASUM kept; lsf_append must refuse
                                    it, lsf_util -repair must restore it bit-exactly, and one
                                    pixel edited then kept with -adopt must change nothing else.
                                    Raises if the test fails.
Stage stats  <p>_noise.csv          one row per file and channel; the .lsf.fits files are loaded
                                    directly (this also tests that Siril opens one as an image).

--export-fits additionally writes <name>.fit next to every LSF product with `lsf_export` (<name>.fit.fz
when FITS compression is enabled in Siril's preferences).

The stack of set 1 is the same frames with the same settings whichever way it is stacked, so it
is computed once (lsf_set_1) and copied to lsf_cum_1 and lsf_1 (siril_set_1 -> siril_cum_1).

NOISE
-----
For every product the script loads it in Siril, runs `stat [-cfa] [main]` and reads the per-channel
statistics back through the Python API. `noise_est_adu` is the value Tools > Image Analysis > Noise
Estimation (the bgnoise command) logs for each channel, and `noise_est_mean_adu` their mean: Siril's
background noise (FnNoise1: the RMS of the differences between neighbouring pixels in each row,
5-sigma clipped, median over the rows, divided by sqrt(2)), x65535 for float images. It comes from
a plain `stat`, run first when --cfa or --stat-main is set. The raw statistic is the `bgnoise`
column. The script also measures an independent per-pixel noise (robust MAD of neighbour
differences on a few horizontal bands; --no-own-noise to skip).

SSIM
----
Structural similarity (Wang et al. 2004: 11x11 Gaussian window, sigma 1.5, K1 0.01, K2 0.03) and
MS-SSIM (Wang, Simoncelli & Bovik 2003: 5 scales, 2x2 averaging between them), pure numpy, one
channel and one strip of rows at a time. Only windows whose pixels all hold data in both images
count (0 or NaN is no data). SSIM's constants assume a known data range, which linear stacks do not
fill, so both images of every pair go through the same transform, set once per channel from
siril_cum_1234 (else the first of the four found): --ssim-domain stretch (default) applies a
Siril-style autostretch (shadows at median - 2.8 x normalized MAD, median to 0.25) and uses range 1;
--ssim-domain linear keeps the data and uses the 0.1-99.9 percentile range. The images share their
frames, and so their noise, so scores sit close to 1: read 1 - SSIM, the 0.1 percentile and the
worst 64x64 blocks (file row and column of the block centre, 1-based) more than the mean.

CLEAN-UP
--------
Scratch files (the stack-on-stack sequence) live in a private sub-folder and are deleted. On ANY
exception the scratch files are removed, Siril's working directory is put back, and your original
frame selection is restored. Finished products are kept.

USAGE
-----
    pyscript lsf_compare.py
    pyscript lsf_compare.py --seq pp_light --outdir D:/tests/m31 --prefix m31
    pyscript lsf_compare.py --rejection sigma --sigma-low 3 --sigma-high 3 --weighting noise
    pyscript lsf_compare.py --split interleaved --skip iter stats
    pyscript lsf_compare.py --dry-run

Needs a Siril build with the lsf_init / lsf_append / lsf_export commands; no other script files.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import itertools
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import traceback

import numpy as np

try:
    import sirilpy as s
    from sirilpy import LogColor
except ImportError:
    print("This script must be run from within Siril (sirilpy not found).")
    raise

N_SETS = 4
STAT_FIELDS = ("total", "ngoodpix", "mean", "median", "sigma", "avgDev", "mad", "sqrtbwmv",
               "location", "scale", "min", "max", "normValue", "bgnoise")
STAGES = ("base", "siril", "sos", "iter", "check", "flip", "parity", "ssim", "stats")   # run order
SKIPPABLE = STAGES

# Siril's 8 rejection methods -> `stack` / lsf_init token, and Siril's default low/high
REJECTION = {
    "none": ("n", None, None),
    "percentile": ("p", 0.2, 0.1),
    "sigma": ("s", 3.0, 3.0),
    "mad": ("m", 3.0, 3.0),
    "sigmedian": ("a", 3.0, 3.0),
    "winsorized": ("w", 3.0, 3.0),
    "linear": ("l", 5.0, 5.0),
    "gesdt": ("g", 0.3, 0.05),
}
NORMS = ("none", "add", "addscale", "mul", "mulscale")
WEIGHTS = ("none", "noise", "wfwhm", "nbstars")
LSF_EXTS = ("LSF_ACC", "LSF_COUNT", "LSF_SETTINGS", "LSF_LEDGER")
SSIM_K1, SSIM_K2 = 0.01, 0.03
SSIM_SIGMA, SSIM_RADIUS = 1.5, 5                                # 11x11 Gaussian window
MSSSIM_WEIGHTS = (0.0448, 0.2856, 0.3001, 0.2363, 0.1333)      # scales 1 (full size) .. 5
SSIM_STRIP = 256                                                # map rows computed at a time
SSIM_BLOCK = 64                                                 # block size of the worst-area report


# --------------------------------------------------------------------------
# Pure helpers (no Siril needed)
# --------------------------------------------------------------------------

def split_sets(indices: list, n_sets: int = N_SETS, mode: str = "contiguous") -> list:
    """Split the selected frame numbers into n_sets sets whose sizes differ by at most one.
    contiguous: consecutive blocks, in sequence order (like successive sessions); the first
    (len % n_sets) sets get the extra frame. interleaved: frame k goes to set k % n_sets."""
    idx = sorted(indices)
    if len(idx) < n_sets:
        raise RuntimeError(f"Need at least {n_sets} selected frames to make {n_sets} sets "
                           f"(the sequence has {len(idx)} selected).")
    if mode == "interleaved":
        return [idx[k::n_sets] for k in range(n_sets)]
    base, extra = divmod(len(idx), n_sets)
    out, pos = [], 0
    for k in range(n_sets):
        size = base + (1 if k < extra else 0)
        out.append(idx[pos:pos + size])
        pos += size
    return out


def runs_of(indices: list) -> list:
    """Sorted unique ints -> list of inclusive (first, last) runs of consecutive numbers."""
    runs = []
    for i in sorted(set(indices)):
        if runs and i == runs[-1][1] + 1:
            runs[-1][1] = i
        else:
            runs.append([i, i])
    return [(a, b) for a, b in runs]


def _parse_card(card: str):
    key = card[:8].strip()
    if card[8:10] != "= ":
        return key, None
    val = card[10:].strip()
    if val.startswith("'"):
        val = val[1:].split("'")[0].rstrip()
    else:
        val = val.split("/")[0].strip()
    return key, val


def read_hdus(path: str) -> list:
    """Walk a FITS file HDU by HDU reading only header blocks. Returns a list of dicts
    {cards: {KEY: raw value}, data_offset, data_size (unpadded bytes)}. Stdlib only."""
    hdus = []
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        pos = 0
        while pos + 2880 <= size:
            f.seek(pos)
            cards, done = {}, False
            while not done:
                block = f.read(2880)
                if len(block) < 2880:
                    raise RuntimeError(f"{os.path.basename(path)}: truncated header at byte {pos}")
                pos += 2880
                for i in range(0, 2880, 80):
                    key, val = _parse_card(block[i:i + 80].decode("ascii", "replace"))
                    if key == "END":
                        done = True
                        break
                    if val is not None and key not in cards:
                        cards[key] = val
            naxis = int(cards.get("NAXIS", 0))
            n = 0
            if naxis:
                n = 1
                for i in range(1, naxis + 1):
                    n *= int(cards.get(f"NAXIS{i}", 0))
            bits = abs(int(cards.get("BITPIX", 8)))
            data = bits // 8 * int(cards.get("GCOUNT", 1)) * (int(cards.get("PCOUNT", 0)) + n)
            if naxis == 0:
                data = 0
            hdus.append({"cards": cards, "data_offset": pos, "data_size": data})
            pos += (data + 2879) // 2880 * 2880
    return hdus


def read_fits_header(path: str) -> dict:
    """Primary-header cards of a FITS file as {KEY: raw value string}."""
    try:
        return read_hdus(path)[0]["cards"]
    except (IndexError, OSError, RuntimeError, ValueError):
        return {}


def header_int(path: str, key: str):
    try:
        return int(float(read_fits_header(path)[key]))
    except (KeyError, ValueError):
        return None


def _q(path: str) -> str:
    return f'"{path}"' if " " in str(path) else str(path)


def _g(x) -> str:
    return f"{float(x):g}"


def plain_path(path: str) -> str:
    """Where the plain-FITS export of an LSF file goes: name.lsf.fits -> name.fit."""
    return (path[:-len(".lsf.fits")] if path.endswith(".lsf.fits") else path) + ".fit"


def saved_path(path: str) -> str:
    """Siril saves a plain image (lsf_export, stack) as path.fz when FITS compression is enabled in
    its preferences: the one of the two that exists, path when neither does."""
    return path + ".fz" if not os.path.exists(path) and os.path.exists(path + ".fz") else path


def remove_saved(path: str) -> None:
    """Removes path and path.fz, so that saved_path cannot find a stale file of an earlier run."""
    for f in (path, path + ".fz"):
        if os.path.exists(f):
            os.remove(f)


def cum_tag(k: int) -> str:
    """Name part of the running stack over sets 1..k: 1, 12, 123, 1234."""
    return "".join(str(x) for x in range(1, k + 1))


def own_noise_from_bands(bands: list, channels: int) -> list:
    """Per-channel noise sigma from horizontal-neighbour differences, pooled over `bands`
    (each an array (channels, h, w) or (h, w)): 1.4826 * MAD(d) / sqrt(2). Differencing removes
    smooth sky and gradients, the MAD ignores stars; zeros (no data) are skipped. None for a
    channel with too little valid data."""
    pools = [[] for _ in range(channels)]
    for b in bands:
        a = np.asarray(b)
        if a.ndim == 2:
            a = a[None]
        if a.shape[0] != channels and a.shape[-1] == channels:
            a = np.moveaxis(a, -1, 0)
        for c in range(min(channels, a.shape[0])):
            x = a[c].astype(np.float64)
            ok = np.isfinite(x) & (x != 0)
            d = (x[:, 1:] - x[:, :-1])[ok[:, 1:] & ok[:, :-1]]
            if d.size:
                pools[c].append(d)
    out = []
    for p in pools:
        d = np.concatenate(p) if p else np.empty(0)
        if d.size < 1000:
            out.append(None)
            continue
        med = float(np.median(d))
        sig = 1.4826 * float(np.median(np.abs(d - med))) / math.sqrt(2.0)
        out.append(sig if (math.isfinite(sig) and sig > 0) else None)
    return out


def as_channels(a) -> np.ndarray:
    """Pixel data as float64 (channels, h, w)."""
    a = np.asarray(a, dtype=np.float64)
    return a[None] if a.ndim == 2 else a


def channel_diffs(x, y) -> list:
    """Per-channel (max abs diff, max rel diff) of x against reference y. The relative
    difference is taken where |y| > 1e-6 (avoids blowing up on empty borders)."""
    x, y = as_channels(x), as_channels(y)
    if x.shape != y.shape:
        raise RuntimeError(f"shape mismatch {x.shape} vs {y.shape}")
    out = []
    for c in range(x.shape[0]):
        d = np.abs(x[c] - y[c])
        ref = np.abs(y[c])
        m = ref > 1e-6
        out.append((float(np.nanmax(d)) if d.size else 0.0,
                    float(np.nanmax(d[m] / ref[m])) if m.any() else 0.0))
    return out


_FITS_DTYPES = {8: ">u1", 16: ">i2", 32: ">i4", 64: ">i8", -32: ">f4", -64: ">f8"}


def fits_image(path: str, extname=None):
    """Image HDU (the primary, or the one named `extname`) read straight from the file as float64
    (channels, h, w) in file row order, with BSCALE/BZERO applied. None when the HDU is missing or
    tile-compressed. Reading every image the same way keeps LSF and Siril files in one orientation."""
    for h in read_hdus(path):
        c = h["cards"]
        if extname is not None and c.get("EXTNAME") != extname:
            continue
        if c.get("ZIMAGE") == "T" or int(c.get("NAXIS", 0)) < 2:
            return None
        dims = [int(c[f"NAXIS{i}"]) for i in range(int(c["NAXIS"]), 0, -1)]
        with open(path, "rb") as f:
            f.seek(h["data_offset"])
            a = np.frombuffer(f.read(h["data_size"]), dtype=_FITS_DTYPES[int(c["BITPIX"])]).reshape(dims)
        a = a.astype(np.float64) * float(c.get("BSCALE", 1)) + float(c.get("BZERO", 0))
        return a[None] if a.ndim == 2 else a
    return None


def noise_estimation(img) -> tuple:
    """Siril's "Noise Estimation" (Tools > Image Analysis, the bgnoise command), ported from
    FnNoise1_float in src/algos/quantize.c: in each row, the first differences between consecutive
    pixels that are not 0 or NaN; their RMS (population sigma), 5-sigma clipped up to 3 times; the
    median over rows, times 1/sqrt(2). Returns (noise, per-row sigma with NaN for skipped rows), in
    the units of `img` (Siril logs it x65535 for float images)."""
    img = np.asarray(img, dtype=np.float32)
    rows = np.full(img.shape[0], np.nan)
    for j, r in enumerate(img):
        v = r[(r != 0) & ~np.isnan(r)]
        if v.size < 3:
            continue
        d = v[:-1] - v[1:]
        dd = d.astype(np.float64)
        mean, sd = dd.mean(), dd.std()
        if sd > 0:
            for _ in range(3):
                keep = np.abs(d - np.float32(mean)) < 5.0 * sd
                if keep.all():
                    break
                d = d[keep]
                dd = d.astype(np.float64)
                mean, sd = dd.mean(), dd.std()
        rows[j] = sd
    ok = rows[~np.isnan(rows)]
    return (0.70710678 * float(np.median(ok)) if ok.size else 0.0), rows


def write_fits_float32(path: str, data, history=()) -> None:
    """Minimal float32 FITS writer (primary HDU only), stdlib + numpy, for diagnostic images."""
    a = np.asarray(data, dtype=">f4")
    if a.ndim == 3 and a.shape[0] == 1:
        a = a[0]
    cards = ["SIMPLE  =                    T", "BITPIX  =                  -32", f"NAXIS   = {a.ndim:20d}"]
    cards += [f"NAXIS{i + 1:<3d}= {n:20d}" for i, n in enumerate(reversed(a.shape))]
    cards += [f"HISTORY {h}"[:80] for h in history] + ["END"]
    head = "".join(c.ljust(80) for c in cards).encode("ascii")
    with open(path, "wb") as f:
        f.write(head + b" " * (-len(head) % 2880))
        raw = a.tobytes()
        f.write(raw + b"\0" * (-len(raw) % 2880))


_TFORM_BYTES = {"L": 1, "X": 1, "B": 1, "I": 2, "J": 4, "K": 8, "A": 1, "E": 4, "D": 8, "C": 8, "M": 16,
                "P": 8, "Q": 16}


def primary_layout(path: str) -> tuple:
    """(data offset, shape (channels, h, w)) of the primary image."""
    c = read_hdus(path)[0]
    cards = c["cards"]
    n = int(cards["NAXIS"])
    shape = (int(cards["NAXIS3"]) if n == 3 else 1, int(cards["NAXIS2"]), int(cards["NAXIS1"]))
    return c["data_offset"], shape


def primary_bytes(path: str) -> bytes:
    off, shape = primary_layout(path)
    with open(path, "rb") as f:
        f.seek(off)
        return f.read(int(np.prod(shape)) * 4)


def edit_primary(path: str, fn) -> None:
    """Rewrites the float32 primary data in place with fn(array (channels, h, w)) -> array, leaving
    every header untouched (so DATASUM and CHECKSUM keep their old values), as an outside tool that
    keeps the extensions would."""
    off, shape = primary_layout(path)
    with open(path, "r+b") as f:
        f.seek(off)
        a = np.frombuffer(f.read(int(np.prod(shape)) * 4), dtype=">f4").reshape(shape)
        f.seek(off)
        f.write(np.ascontiguousarray(fn(a.copy()), dtype=">f4").tobytes())


def parity_breakdown(x, y, count=None, stackcnt=None) -> list:
    """Per-channel parity of an LSF stack x against a Siril stack y (both (channels, h, w), same
    orientation). Splits the pixels into full coverage (LSF_COUNT == STACKCNT), partial coverage,
    and pixels where exactly one image is 0 (no data in one of them), and compares Siril's noise
    estimate on both, row by row."""
    x, y = as_channels(x), as_channels(y)
    if x.shape != y.shape:
        raise RuntimeError(f"shape mismatch {x.shape} vs {y.shape}")
    out = []
    for c in range(x.shape[0]):
        d = np.abs(x[c] - y[c])
        one_zero = (x[c] == 0) ^ (y[c] == 0)
        both = ~one_zero
        r = {"n_diff": int((d > 0).sum()), "max_abs": float(d.max()) if d.size else 0.0,
             "max_abs_both_data": float(d[both].max()) if both.any() else 0.0,
             "n_one_zero": int(one_zero.sum())}
        if count is not None and stackcnt:
            n = count[c] if count.shape[0] > c else count[0]
            full = n == stackcnt
            r["n_diff_full"] = int((d[full] > 0).sum())
            r["n_diff_partial"] = int((d[~full] > 0).sum())
            r["max_abs_full"] = float(d[full].max()) if full.any() else 0.0
            r["max_abs_partial"] = float(d[~full].max()) if (~full).any() else 0.0
        nx, rx = noise_estimation(x[c])
        ny, ry = noise_estimation(y[c])
        rel = np.abs(rx - ry) / np.where(ry > 0, ry, np.nan)
        r["noise_lsf"], r["noise_siril"] = nx, ny
        r["noise_ratio"] = nx / ny if ny else float("nan")
        r["rows_noise_differ"] = int(np.nansum(rel > 1e-4))
        worst = np.argsort(np.nan_to_num(rel, nan=-1.0))[::-1][:5]
        r["worst_rows"] = [(int(j), float(rel[j])) for j in worst if np.isfinite(rel[j]) and rel[j] > 1e-4]
        out.append(r)
    return out


def fits_memmap(path: str):
    """The primary image as (read-only memory map (channels, h, w) in file row order, BSCALE, BZERO),
    so that one channel at a time is read; None when it is tile-compressed or not an image."""
    h = read_hdus(path)[0]
    c = h["cards"]
    if c.get("ZIMAGE") == "T" or int(c.get("NAXIS", 0)) < 2:
        return None
    dims = tuple(int(c[f"NAXIS{i}"]) for i in range(int(c["NAXIS"]), 0, -1))
    mm = np.memmap(path, dtype=_FITS_DTYPES[int(c["BITPIX"])], mode="r", offset=h["data_offset"], shape=dims)
    return (mm[None] if mm.ndim == 2 else mm), float(c.get("BSCALE", 1)), float(c.get("BZERO", 0))


def channel_f32(img, c: int) -> np.ndarray:
    """Channel c of an image from fits_memmap (or (array, 1, 0)) as a new native float32 array."""
    a, bscale, bzero = img
    x = np.array(a[c], dtype=np.float32)
    if bscale != 1:
        x *= np.float32(bscale)
    if bzero:
        x += np.float32(bzero)
    return x


def mtf(x, m):
    """Midtones transfer function: 0 -> 0, m -> 0.5, 1 -> 1."""
    return (m - 1) * x / ((2 * m - 1) * x - m)


def ssim_reference(values, domain: str) -> tuple:
    """(shadows, midtones, data range L) for one channel, from the valid pixels of the reference
    image. stretch: Siril-style autostretch, shadows at median - 2.8 x normalized MAD and the median
    taken to 0.25, L = 1. linear: no transform (None, None) and L = 99.9th - 0.1th percentile."""
    v = np.asarray(values[::max(1, values.size // 2_000_000)], dtype=np.float64)
    if not v.size:
        raise RuntimeError("the reference channel has no data")
    if domain == "linear":
        lo, hi = np.percentile(v, (0.1, 99.9))
        return None, None, float(hi - lo) if hi > lo else 1.0
    med = float(np.median(v))
    mad = 1.4826 * float(np.median(np.abs(v - med)))
    c0 = min(max(med - 2.8 * mad, 0.0), med)
    b = (med - c0) / (1 - c0) if c0 < 1 else 0.0
    return c0, (float(mtf(b, 0.25)) if 0 < b < 1 else 0.5), 1.0


def apply_stretch(x: np.ndarray, c0, m) -> np.ndarray:
    """ssim_reference's transform on float32 x (changed in place for the shadows step)."""
    if c0 is None:
        return x
    x -= np.float32(c0)
    x *= np.float32(1 / (1 - c0))
    np.clip(x, 0, 1, out=x)
    return mtf(x, np.float32(m)).astype(np.float32, copy=False)


def _filter_valid(a: np.ndarray, taps: np.ndarray) -> np.ndarray:
    """Separable correlation of the last two axes with `taps`, 'valid' mode: each shrinks by
    len(taps) - 1."""
    n = len(taps)
    w = a.shape[-1] - n + 1
    h = a[..., :w] * taps[0]
    tmp = np.empty_like(h)
    for k in range(1, n):
        np.multiply(a[..., k:k + w], taps[k], out=tmp)
        h += tmp
    rows = h.shape[-2] - n + 1
    v = h[..., :rows, :] * taps[0]
    tmp = tmp[..., :rows, :]
    for k in range(1, n):
        np.multiply(h[..., k:k + rows, :], taps[k], out=tmp)
        v += tmp
    return v


def ssim_scale(x, y, valid, c1: float, c2: float, keep_map: bool = False):
    """SSIM of two float32 images (h, w) at one scale, strip by strip: (mean SSIM, mean
    contrast-structure term, number of windows, SSIM map or None). Only windows that lie wholly on
    `valid` pixels count; the map is (h - 10, w - 10), window centres from pixel (5, 5), NaN at the
    other windows. None when the image is smaller than the window."""
    n = 2 * SSIM_RADIUS + 1
    hh, ww = x.shape[0] - n + 1, x.shape[1] - n + 1
    if hh < 1 or ww < 1:
        return None
    k = np.arange(-SSIM_RADIUS, SSIM_RADIUS + 1)
    taps = np.exp(-k * k / (2 * SSIM_SIGMA ** 2))
    taps = (taps / taps.sum()).astype(np.float32)
    box = np.ones(n, np.float32)
    smap = np.full((hh, ww), np.nan, np.float32) if keep_map else None
    s_sum = cs_sum = 0.0
    count = 0
    for r0 in range(0, hh, SSIM_STRIP):
        r1 = min(hh, r0 + SSIM_STRIP)
        xs, ys, vs = x[r0:r1 + n - 1], y[r0:r1 + n - 1], valid[r0:r1 + n - 1]
        mx, my, xx, yy, xy = _filter_valid(np.stack((xs, ys, xs * xs, ys * ys, xs * ys)), taps)
        ok = (np.ones((r1 - r0, ww), bool) if vs.all()
              else _filter_valid(vs.astype(np.float32), box) > n * n - 0.5)
        mxy = mx * my
        mx *= mx
        my *= my
        cs = (2 * (xy - mxy) + c2) / ((xx - mx) + (yy - my) + c2)
        s = (2 * mxy + c1) / (mx + my + c1) * cs
        count += int(ok.sum())
        s_sum += float(s[ok].sum(dtype=np.float64))
        cs_sum += float(cs[ok].sum(dtype=np.float64))
        if smap is not None:
            smap[r0:r1][ok] = s[ok]
    if not count:
        return float("nan"), float("nan"), 0, smap
    return s_sum / count, cs_sum / count, count, smap


def _halve(x, valid) -> tuple:
    """2x2 average (MS-SSIM's downsampling); a pixel keeps data only if all four had it."""
    h, w = x.shape[0] // 2 * 2, x.shape[1] // 2 * 2
    xd = (x[0:h:2, 0:w:2] + x[1:h:2, 0:w:2] + x[0:h:2, 1:w:2] + x[1:h:2, 1:w:2]) * np.float32(0.25)
    vd = valid[0:h:2, 0:w:2] & valid[1:h:2, 0:w:2] & valid[0:h:2, 1:w:2] & valid[1:h:2, 1:w:2]
    return xd, vd


def ms_ssim(x, y, valid, c1: float, c2: float, first) -> float:
    """MS-SSIM: contrast-structure at scales 1-4 and full SSIM at scale 5, weighted by MSSSIM_WEIGHTS
    (terms below 0 count as 0). `first` is ssim_scale's result at full size. NaN when the image is
    too small for 5 scales or a scale has no valid window."""
    if first is None or not first[2]:
        return float("nan")
    terms = [first[1]]
    last = len(MSSSIM_WEIGHTS) - 1
    for j in range(1, last + 1):
        (x, valid), (y, _) = _halve(x, valid), _halve(y, valid)
        r = ssim_scale(x, y, valid, c1, c2)
        if r is None or not r[2]:
            return float("nan")
        terms.append(r[0] if j == last else r[1])
    return float(np.prod([max(t, 0.0) ** wt for t, wt in zip(terms, MSSSIM_WEIGHTS)]))


def worst_blocks(smap, size: int = SSIM_BLOCK, k: int = 3) -> list:
    """The k SSIM_BLOCK-square blocks of the SSIM map with the lowest mean (blocks with at least half
    their windows valid), as (file row, file column, mean SSIM), centre pixel, 1-based."""
    hh, ww = smap.shape[0] // size * size, smap.shape[1] // size * size
    if not hh or not ww:
        return []
    b = smap[:hh, :ww].reshape(hh // size, size, ww // size, size)
    ok = np.isfinite(b)
    cnt = ok.sum(axis=(1, 3))
    mean = np.where(cnt >= size * size // 2,
                    np.where(ok, b, 0).sum(axis=(1, 3), dtype=np.float64) / np.maximum(cnt, 1), np.inf)
    out = []
    for i in np.argsort(mean, axis=None)[:k]:
        bi, bj = divmod(int(i), mean.shape[1])
        if np.isfinite(mean[bi, bj]):
            out.append((bi * size + size // 2 + SSIM_RADIUS + 1, bj * size + size // 2 + SSIM_RADIUS + 1,
                        float(mean[bi, bj])))
    return out


def ssim_channel(x, y, transform: tuple, keep_map: bool = False) -> dict:
    """SSIM, its 0.1 percentile and minimum over the windows, MS-SSIM and the worst blocks of one
    channel pair (float32 (h, w), consumed). transform = ssim_reference(...). With keep_map, also
    "dissim": 1 - SSIM at the window centres, image size, 0 where no window counts."""
    c0, m, L = transform
    valid = np.isfinite(x) & np.isfinite(y) & (x != 0) & (y != 0)
    x, y = apply_stretch(x, c0, m), apply_stretch(y, c0, m)
    x[~valid] = 0
    y[~valid] = 0
    c1, c2 = (SSIM_K1 * L) ** 2, (SSIM_K2 * L) ** 2
    first = ssim_scale(x, y, valid, c1, c2, keep_map=True)
    if first is None:
        raise RuntimeError(f"image {x.shape} is smaller than the SSIM window")
    smap = first[3]
    v = smap[np.isfinite(smap)]
    r = {"ssim": first[0], "windows": first[2], "data_range": L,
         "ssim_p0_1": float(np.percentile(v, 0.1)) if v.size else float("nan"),
         "ssim_min": float(v.min()) if v.size else float("nan")}
    del v
    r["ms_ssim"] = ms_ssim(x, y, valid, c1, c2, first)
    r["worst"] = worst_blocks(smap)
    if keep_map:
        d = np.zeros(x.shape, np.float32)
        d[SSIM_RADIUS:SSIM_RADIUS + smap.shape[0], SSIM_RADIUS:SSIM_RADIUS + smap.shape[1]] = \
            np.nan_to_num(1 - smap, nan=0.0)
        r["dissim"] = d
    return r


def check_lsf_file(path: str) -> tuple:
    """v0.4 structural checks of one LSF file. Returns (results, info) where results is a list of
    (check, ok, detail) and info = {uid, gen}."""
    res = []

    def chk(name, ok, detail=""):
        res.append((name, bool(ok), detail))

    hdus = read_hdus(path)
    prim = hdus[0]["cards"]
    chk("primary_lsf", prim.get("LSF") == "T", f"LSF={prim.get('LSF')}")
    chk("lsfver", prim.get("LSFVER") == "0.4", f"LSFVER={prim.get('LSFVER')}")
    chk("primary_float32", prim.get("BITPIX") == "-32", f"BITPIX={prim.get('BITPIX')}")
    names = [h["cards"].get("EXTNAME", "") for h in hdus[1:]]
    for e in LSF_EXTS:
        chk(f"has_{e}", e in names, ",".join(names))
    chk("ledger_last", bool(names) and names[-1] == "LSF_LEDGER", f"last={names[-1] if names else None}")
    by = {h["cards"].get("EXTNAME", ""): h["cards"] for h in hdus[1:]}
    p = int(prim.get("NAXIS3", 1)) if int(prim.get("NAXIS", 2)) >= 3 else 1
    acc = by.get("LSF_ACC", {})
    n3 = int(acc.get("NAXIS3", 0)) if acc else 0
    # 5*P when the file has an exposure map (EXP planes after W2)
    chk("acc_naxis3_4P", n3 in (4 * p, 5 * p), f"NAXIS3={n3} P={p}")
    planes = [acc.get(f"LSFPL{i}") for i in range(1, n3 + 1)]
    chk("acc_planes_named", all(planes), ",".join(str(x) for x in planes))
    if n3 == 5 * p:
        want = ["EXP"] if p == 1 else ["EXP_R", "EXP_G", "EXP_B"]
        chk("acc_exp_planes", all(w in [str(x).strip() for x in planes] for w in want), ",".join(str(x) for x in planes))
    uid, gen = prim.get("LSFUID"), prim.get("LSFGEN")
    chk("primary_uid_gen", uid is not None and gen is not None, f"LSFUID={uid} LSFGEN={gen}")
    for e in LSF_EXTS:
        c = by.get(e)
        if c is not None:
            chk(f"{e}_uid_gen", c.get("LSFUID") == uid and c.get("LSFGEN") == gen,
                f"LSFUID={c.get('LSFUID')} LSFGEN={c.get('LSFGEN')}")
    led = by.get("LSF_LEDGER", {})
    stackcnt, rows = prim.get("STACKCNT"), led.get("NAXIS2")
    chk("stackcnt_eq_ledger", stackcnt is not None and rows is not None and int(float(stackcnt)) == int(rows),
        f"STACKCNT={stackcnt} ledger rows={rows}")
    with open(path, "rb") as f:
        f.seek(hdus[0]["data_offset"])
        h = hashlib.sha256()
        left = hdus[0]["data_size"]
        while left:
            buf = f.read(min(left, 1 << 22))
            if not buf:
                break
            h.update(buf)
            left -= len(buf)
    want = (prim.get("LSFHASH") or "").lower()
    chk("lsfhash", h.hexdigest() == want, f"stored={want[:16]}... computed={h.hexdigest()[:16]}...")
    try:
        gen_i = int(float(gen)) if gen is not None else None
    except ValueError:
        gen_i = None
    return res, {"uid": uid, "gen": gen_i}


def var_string_bad_rows(path: str, hdu_index: int, col: int) -> list:
    """Rows (1-based) whose variable-length string (TFORM nPA / nQA) in column `col` of HDU
    `hdu_index` (0-based) holds bytes outside printable ASCII, read straight from the heap."""
    hdu = read_hdus(path)[hdu_index]
    c = hdu["cards"]
    width, nrows = int(c["NAXIS1"]), int(c["NAXIS2"])
    offs, pos, kind = {}, 0, None
    for i in range(1, int(c["TFIELDS"]) + 1):
        m = re.match(r"\s*(\d*)([A-Z])(\w?)", c[f"TFORM{i}"])
        rep, code = int(m.group(1) or 1), m.group(2)
        offs[i] = pos
        if i == col:
            kind = code if code in "PQ" and m.group(3) == "A" else None
        pos += _TFORM_BYTES[code] * (1 if code in "PQ" else rep)
    if kind is None:
        raise ValueError(f"column {col} is not a variable-length string column")
    heap = hdu["data_offset"] + int(c.get("THEAP", width * nrows))
    bad = []
    with open(path, "rb") as f:
        for r in range(nrows):
            f.seek(hdu["data_offset"] + r * width + offs[col])
            n, off = struct.unpack(">ii" if kind == "P" else ">qq", f.read(8 if kind == "P" else 16))
            f.seek(heap + off)
            if any(b < 32 or b > 126 for b in f.read(n)):
                bad.append(r + 1)
    return bad


def run_fitsverify(path: str):
    """(ok, detail) from fitsverify, or None if it is not on PATH.

    CFITSIO (seen in 4.7.0) truncates a variable-length string when reading it across two FITS block
    boundaries, and fitsverify then scans stale bytes of its reused buffer and reports "String in row
    #r, and column #c contains non-ASCII text". Such errors are re-checked byte for byte from the heap
    and ignored when the stored string is clean; any other fitsverify error is a failure."""
    exe = shutil.which("fitsverify")
    if not exe:
        return None
    r = subprocess.run([exe, path], capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip()
    m = re.search(r"found\s+(\d+)\s+warning\(s\) and\s+(\d+)\s+error", out)
    errors = int(m.group(2)) if m else (0 if r.returncode == 0 else 1)
    if errors == 0:
        return True, "verification OK"
    hdu, other, suspect = 0, [], {}
    for line in out.splitlines():
        h = re.search(r"=+\s*HDU\s+(\d+):", line)
        if h:
            hdu = int(h.group(1)) - 1                     # fitsverify numbers HDUs from 1
        elif "*** Error:" in line:
            v = re.search(r"String in row #(\d+), and column #(\d+) contains non-ASCII text", line)
            if v:
                suspect.setdefault((hdu, int(v.group(2))), []).append(int(v.group(1)))
            else:
                other.append(line.split("*** Error:", 1)[1].strip())
    if "give up" in out:
        other.append("fitsverify gave up (too many errors)")
    n_suspect = sum(len(v) for v in suspect.values())
    if n_suspect + len(other) != errors:
        other.append(f"{errors} errors reported, {n_suspect + len(other)} parsed")
    for (h, col), rows in suspect.items():
        try:
            bad = var_string_bad_rows(path, h, col)
        except Exception as e:  # noqa: BLE001
            other.append(f"HDU {h + 1} col {col}: cannot re-check ({e})")
            continue
        if bad:
            other.append(f"HDU {h + 1} col {col}: non-ASCII bytes in rows {bad[:10]}")
    if other:
        return False, f"{errors} error(s): " + "; ".join(other[:5])
    return True, (f"verification OK apart from {n_suspect} CFITSIO long-string read false positive(s) "
                  f"(stored strings checked byte for byte)")


def run_cmd(siril, *parts) -> None:
    """siril.cmd that always raises RuntimeError on failure (whether it returns False or raises)."""
    text = " ".join(str(p) for p in parts)
    try:
        ok = siril.cmd(*[str(p) for p in parts])
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"Siril command failed: {text} ({e})") from e
    if ok is False:
        raise RuntimeError(f"Siril command failed: {text}")


# --------------------------------------------------------------------------
# Arguments
# --------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="lsf_compare.py",
        description="End-to-end test of Siril's lsf_init / lsf_append / lsf_export against Siril's stack.")
    p.add_argument("--seq", default=None, help="Sequence to use (default: the sequence loaded in Siril).")
    p.add_argument("--outdir", default=None, help="Output folder (default: <Siril working dir>/lsf_compare).")
    p.add_argument("--prefix", default="cmp", help="File name prefix (default: cmp).")
    p.add_argument("--split", choices=["contiguous", "interleaved"], default="contiguous",
                   help="How the selected frames are divided into the 4 sets (default: contiguous blocks).")
    p.add_argument("--overwrite", action="store_true", help="Replace existing products instead of refusing to run.")
    p.add_argument("--dry-run", action="store_true", help="Print the plan (sets, files, commands) and stop.")
    p.add_argument("--skip", nargs="+", choices=SKIPPABLE, default=[],
                   help="Leave out stages: base (lsf_init per-set + running stacks), siril (native stack), "
                        "sos (stack on stack), iter (lsf_append chain), check (v0.4 structure), "
                        "flip (repair/adopt test), parity (LSF vs Siril), ssim (SSIM of the full stacks), "
                        "stats (noise export).")
    p.add_argument("--export-fits", action="store_true",
                   help="Also write <name>.fit next to every LSF file with lsf_export.")
    g = p.add_argument_group("stacking settings (shared by lsf_init and stack)")
    g.add_argument("--rejection", choices=list(REJECTION), default="winsorized")
    g.add_argument("--sigma-low", type=float, default=None, help="Default depends on the method (Siril's).")
    g.add_argument("--sigma-high", type=float, default=None, help="Default depends on the method (Siril's).")
    g.add_argument("--weighting", choices=list(WEIGHTS), default="none")
    g.add_argument("--norm", choices=list(NORMS), default="addscale")
    g.add_argument("--fastnorm", action="store_true", help="Pass -fastnorm.")
    g.add_argument("--rgb-equal", action=argparse.BooleanOptionalAction, default=True,
                   help="Colour: pass -rgb_equal. Default on; ignored for mono sequences.")
    g.add_argument("--sos-rejection", choices=list(REJECTION), default=None,
                   help="Rejection for the stack on stack only (default: same as --rejection).")
    g.add_argument("--sos-sigma-low", type=float, default=None)
    g.add_argument("--sos-sigma-high", type=float, default=None)
    g = p.add_argument_group("noise export")
    g.add_argument("--cfa", action="store_true", help="Add -cfa to the stat command (undebayered data).")
    g.add_argument("--stat-main", action="store_true",
                   help="Use `stat main` (faster, but without Siril's noise estimate).")
    g.add_argument("--no-own-noise", dest="own_noise", action="store_false", default=True,
                   help="Skip the independent neighbour-difference noise estimate.")
    g = p.add_argument_group("SSIM")
    g.add_argument("--ssim-domain", choices=["stretch", "linear"], default="stretch",
                   help="Compare autostretched images (range 1, default) or the linear data (range = "
                        "0.1-99.9 percentile of siril_cum_1234).")
    g.add_argument("--ssim-maps", action="store_true", help="Also write a 1 - SSIM map for every pair.")
    return p


# --------------------------------------------------------------------------
# The comparison run
# --------------------------------------------------------------------------

class Comparison:
    def __init__(self, siril, args):
        self.siril, self.a = siril, args
        self.log = siril.log
        self.seqname = None
        self.nframes = 0
        self.original_selected = []
        self.sets = []
        self.wd = None
        self.outdir = None
        self.workdir = None
        self.products = []              # (label, kind, path, expected frame count or None)
        self.done = {}                  # label -> path of products made in THIS run
        self._own_noise_ok = True
        self.selection_touched = False  # only restore the user's selection if we ever changed it
        self.t0 = time.perf_counter()

    # ---- paths -------------------------------------------------------------
    def p(self, suffix: str) -> str:
        return os.path.join(self.outdir, f"{self.a.prefix}_{suffix}")

    def planned_files(self) -> list:
        a = self.a
        tags = [cum_tag(k) for k in range(1, N_SETS + 1)]
        out = []
        if "base" not in a.skip:
            out += [self.p(f"lsf_set_{k}.lsf.fits") for k in range(1, N_SETS + 1)]
            out += [self.p(f"lsf_cum_{t}.lsf.fits") for t in tags]
        if "siril" not in a.skip:
            out += [self.p(f"siril_set_{k}.fit") for k in range(1, N_SETS + 1)]
            out += [self.p(f"siril_cum_{t}.fit") for t in tags]
        if "sos" not in a.skip:
            out.append(self.p("lsf_sos.lsf.fits"))
        if "iter" not in a.skip:
            out += [self.p(f"lsf_{k}.lsf.fits") for k in range(1, N_SETS + 1)]
        if a.export_fits:
            out += [plain_path(f) for f in out if f.endswith(".lsf.fits")]
        for st, name in (("check", "check.csv"), ("flip", "flip.csv"), ("parity", "parity.csv"),
                         ("ssim", "ssim.csv"), ("stats", "noise.csv")):
            if st not in a.skip:
                out.append(self.p(name))
        return out

    # ---- settings -> command arguments ----------------------------------------
    def rej_args(self, rejection=None, low=None, high=None) -> list:
        a = self.a
        rej = rejection or a.rejection
        if rejection is None:
            low, high = a.sigma_low, a.sigma_high
        tok, dlo, dhi = REJECTION[rej]
        if rej == "none":
            return ["rej", tok]
        return ["rej", tok, _g(dlo if low is None else low), _g(dhi if high is None else high)]

    def common_opts(self) -> list:
        a = self.a
        out = []
        if a.weighting != "none":
            out.append(f"-weight={a.weighting}")
        out.append("-nonorm" if a.norm == "none" else f"-norm={a.norm}")
        if a.fastnorm and a.norm != "none":
            out.append("-fastnorm")
        if a.rgb_equal:
            out.append("-rgb_equal")
        return out

    def lsf_init_cmd(self, dest: str, seq=None, rejection=None, low=None, high=None) -> list:
        return (["lsf_init", _q(seq or self.seqname), _q(dest)] + self.rej_args(rejection, low, high)
                + self.common_opts() + ["-overwrite"])

    def lsf_append_cmd(self, path: str) -> list:
        return ["lsf_append", _q(path), _q(self.seqname), "-noreg"]

    def siril_stack_cmd(self, dest: str) -> list:
        return (["stack", _q(self.seqname)] + self.rej_args() + self.common_opts()
                + ["-32b", "-filter-included", _q(f"-out={dest}")])

    # ---- sequence / selection ----------------------------------------------
    def current_seq(self):
        try:
            return self.siril.get_seq()
        except Exception:  # noqa: BLE001
            return None

    def ensure_seq(self) -> None:
        """Make sure our sequence is the loaded one (loading single images closes it)."""
        seq = self.current_seq()
        if seq is None or seq.seqname != self.seqname:
            run_cmd(self.siril, "load_seq", _q(self.seqname))

    def apply_selection(self, wanted: list) -> None:
        """Make exactly the frames `wanted` (0-based) the selected ones, using Siril's own
        select / unselect commands (they keep the .seq file in sync), then verify."""
        self.ensure_seq()
        self.selection_touched = True
        n, name = self.nframes, _q(self.seqname)
        wanted = sorted(set(wanted))
        wset = set(wanted)
        sel_runs = runs_of(wanted)
        gap_runs = runs_of([i for i in range(n) if i not in wset])
        if len(sel_runs) <= len(gap_runs):
            run_cmd(self.siril, "unselect", name, 1, n)
            for a, b in sel_runs:
                run_cmd(self.siril, "select", name, a + 1, b + 1)
        else:
            run_cmd(self.siril, "select", name, 1, n)
            for a, b in gap_runs:
                run_cmd(self.siril, "unselect", name, a + 1, b + 1)
        seq = self.current_seq()
        got = [i for i in range(seq.number) if seq.imgparam[i].incl]
        if got != wanted:
            raise RuntimeError(f"Frame selection did not take: wanted {len(wanted)} frames, Siril reports "
                               f"{len(got)} selected.")

    # ---- bookkeeping -----------------------------------------------------------
    def check_count(self, path: str, expected, label: str) -> None:
        got = header_int(path, "STACKCNT")
        if expected is not None and got is not None and got != expected:
            raise RuntimeError(f"{label}: STACKCNT is {got} but {expected} frames were meant to be stacked "
                               f"(the frame selection was not honoured).")

    def register(self, label, kind, path, expected) -> None:
        path = saved_path(path)
        if not os.path.exists(path):
            raise RuntimeError(f"{label}: the command reported success but {os.path.basename(path)} is missing.")
        self.check_count(path, expected, label)
        self.products.append((label, kind, path, expected))
        self.done[label] = path
        if kind.startswith("lsf") and self.a.export_fits:
            self.export_plain(path)

    def export_plain(self, lsf_path: str) -> None:
        """lsf_export of an LSF file to <name>.fit (.fit.fz with Siril's FITS compression). A failed
        export is a warning."""
        dest = plain_path(lsf_path)
        try:
            remove_saved(dest)
            run_cmd(self.siril, "lsf_export", _q(lsf_path), _q(dest))
            self.log(f"  lsf_export -> {os.path.basename(saved_path(dest))}")
        except Exception as e:  # noqa: BLE001
            self.log(f"WARNING: lsf_export of {os.path.basename(lsf_path)} failed ({e}).", LogColor.SALMON)

    def copy_product(self, src_label: str, dest: str, label: str, kind: str, expected, why: str) -> None:
        """A product that is the very same stack as an earlier one: copy instead of recomputing."""
        if os.path.exists(dest):
            os.remove(dest)
        shutil.copy2(self.done[src_label], dest)
        self.register(label, kind, dest, expected)
        self.log(f"{label}: {why}; copied from {src_label} instead of recomputed.")

    def timed(self, cmd: list, label: str, what: str) -> None:
        t = time.perf_counter()
        run_cmd(self.siril, *cmd)
        self.log(f"{label}: {what} in {time.perf_counter() - t:.1f} s", LogColor.GREEN)

    # ---- stages -------------------------------------------------------------
    def stage_base(self, step) -> None:
        """lsf_init per-set stacks, then running stacks over sets 1..k, each fresh."""
        S = self.sets
        for k in range(1, N_SETS + 1):
            step(f"lsf_init of set {k}")
            dest = self.p(f"lsf_set_{k}.lsf.fits")
            self.apply_selection(S[k - 1])
            self.timed(self.lsf_init_cmd(dest), f"lsf_set_{k}", f"{len(S[k - 1])} frames stacked")
            self.register(f"lsf_set_{k}", "lsf_base", dest, len(S[k - 1]))
        self.copy_product("lsf_set_1", self.p("lsf_cum_1.lsf.fits"), "lsf_cum_1", "lsf_base", len(S[0]),
                          "the running stack over set 1 IS the stack of set 1")
        for k in range(2, N_SETS + 1):
            step(f"lsf_init of sets 1-{k}")
            frames = sorted(i for grp in S[:k] for i in grp)
            dest = self.p(f"lsf_cum_{cum_tag(k)}.lsf.fits")
            self.apply_selection(frames)
            self.timed(self.lsf_init_cmd(dest), f"lsf_cum_{cum_tag(k)}", f"{len(frames)} frames stacked")
            self.register(f"lsf_cum_{cum_tag(k)}", "lsf_base", dest, len(frames))

    def stage_siril(self, step) -> None:
        """Native Siril stack of each set and of each running set, same settings."""
        S = self.sets
        for k in range(1, N_SETS + 1):
            step(f"Siril stack of set {k}")
            dest = self.p(f"siril_set_{k}.fit")
            remove_saved(dest)
            self.apply_selection(S[k - 1])
            self.timed(self.siril_stack_cmd(dest), f"siril_set_{k}", f"{len(S[k - 1])} frames stacked")
            self.register(f"siril_set_{k}", "siril", dest, len(S[k - 1]))
        self.copy_product("siril_set_1", self.p("siril_cum_1.fit"), "siril_cum_1", "siril", len(S[0]),
                          "the running stack over set 1 IS the stack of set 1")
        for k in range(2, N_SETS + 1):
            step(f"Siril stack of sets 1-{k}")
            frames = sorted(i for grp in S[:k] for i in grp)
            dest = self.p(f"siril_cum_{cum_tag(k)}.fit")
            remove_saved(dest)
            self.apply_selection(frames)
            self.timed(self.siril_stack_cmd(dest), f"siril_cum_{cum_tag(k)}", f"{len(frames)} frames stacked")
            self.register(f"siril_cum_{cum_tag(k)}", "siril", dest, len(frames))

    def stage_iter(self, step) -> None:
        """lsf_1 = set 1; each next file = copy of the previous one, then lsf_append of the next set."""
        total, prev = 0, None
        for k in range(1, N_SETS + 1):
            step(f"LSF {'initialize from set 1' if k == 1 else f'append set {k}'}")
            frames = self.sets[k - 1]
            total += len(frames)
            dest = self.p(f"lsf_{k}.lsf.fits")
            if k == 1:
                if "lsf_set_1" in self.done:
                    self.copy_product("lsf_set_1", dest, "lsf_1", "lsf_iter", total,
                                      "initializing from set 1 is the stack of set 1")
                else:
                    self.apply_selection(frames)
                    self.timed(self.lsf_init_cmd(dest), "lsf_1", f"{len(frames)} frames initialized")
                    self.register("lsf_1", "lsf_iter", dest, total)
            else:
                if os.path.abspath(prev) != os.path.abspath(dest):   # lsf_append writes in place
                    if os.path.exists(dest):
                        os.remove(dest)
                    shutil.copy2(prev, dest)
                self.apply_selection(frames)
                self.timed(self.lsf_append_cmd(dest), f"lsf_{k}",
                           f"{len(frames)} frames appended ({total} in history)")
                self.register(f"lsf_{k}", "lsf_iter", dest, total)
            prev = dest

    def stage_sos(self, step) -> None:
        """LSF stack on stack: the four per-set LSF stacks (lsf_export'ed) become a scratch sequence,
        stacked with lsf_init. Scratch is removed straight away."""
        a = self.a
        srcs = [self.p(f"lsf_set_{k}.lsf.fits") for k in range(1, N_SETS + 1)]
        missing = [x for x in srcs if not os.path.exists(x)]
        if missing:
            raise RuntimeError("LSF stack on stack needs the per-set LSF stacks; missing: "
                               + ", ".join(os.path.basename(m) for m in missing))
        step("lsf_init of the four per-set LSF stacks")
        sos_dir = tempfile.mkdtemp(prefix="sos_", dir=self.workdir)
        root, cd_done = "lsfsos", False
        try:
            for i, src in enumerate(srcs, 1):
                run_cmd(self.siril, "lsf_export", _q(src), _q(os.path.join(sos_dir, f"{root}{i:05d}.fit")))
            run_cmd(self.siril, "cd", _q(sos_dir))
            cd_done = True
            if not self.siril.create_new_seq(root):
                raise RuntimeError("Siril could not create the stack-on-stack sequence.")
            name = next((c for c in (root, root + "_") if os.path.exists(os.path.join(sos_dir, c + ".seq"))), root)
            run_cmd(self.siril, "load_seq", _q(name))
            dest = self.p("lsf_sos.lsf.fits")
            self.timed(self.lsf_init_cmd(dest, name, a.sos_rejection, a.sos_sigma_low, a.sos_sigma_high),
                       "lsf_sos", f"4 per-set stacks stacked (rejection {a.sos_rejection or a.rejection})")
        finally:
            if cd_done:
                try:
                    run_cmd(self.siril, "cd", _q(self.wd))
                except Exception as e:  # noqa: BLE001
                    self.log(f"Could not restore Siril's working directory: {e}", LogColor.SALMON)
        self.register("lsf_sos", "lsf_base", dest, N_SETS)
        shutil.rmtree(sos_dir, ignore_errors=True)

    def stage_check(self, step) -> None:
        """Pure-Python v0.4 structure checks of every .lsf.fits product; raises on any failure."""
        step("LSF v0.4 structure checks")
        files = [(lbl, self.p(name)) for lbl, kind, name in self.known_outputs()
                 if name.endswith(".lsf.fits") and os.path.exists(self.p(name))]
        rows, failed, gens = [], [], {}
        for label, path in files:
            try:
                res, info = check_lsf_file(path)
            except Exception as e:  # noqa: BLE001
                res, info = [("read", False, str(e))], {"gen": None}
            fv = run_fitsverify(path)
            if fv is not None:
                res.append(("fitsverify", fv[0], fv[1]))
            gens[label] = info.get("gen")
            for name, ok, detail in res:
                rows.append((label, os.path.basename(path), name, "ok" if ok else "FAIL", detail))
                if not ok:
                    failed.append(f"{label}:{name} ({detail})")
        for k in range(2, N_SETS + 1):
            a_, b_ = gens.get(f"lsf_{k - 1}"), gens.get(f"lsf_{k}")
            if f"lsf_{k}" in gens:
                ok = a_ is not None and b_ is not None and b_ == a_ + 1
                rows.append((f"lsf_{k}", "", "lsfgen_increment", "ok" if ok else "FAIL",
                             f"lsf_{k - 1}={a_} lsf_{k}={b_}"))
                if not ok:
                    failed.append(f"lsf_{k}:lsfgen_increment ({a_} -> {b_})")
        if shutil.which("fitsverify") is None:
            self.log("fitsverify not on PATH: skipped.", LogColor.SALMON)
        csv_path = self.p("check.csv")
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["label", "file", "check", "result", "detail"])
            w.writerows(rows)
        self.log(f"Structure checks: {len(rows) - len(failed)} ok, {len(failed)} failed, on {len(files)} "
                 f"files -> {csv_path}", LogColor.RED if failed else LogColor.GREEN)
        if failed:
            raise RuntimeError("LSF v0.4 structure checks failed: " + "; ".join(failed[:10])
                               + (" ..." if len(failed) > 10 else ""))

    def known_exceptions(self) -> list:
        a = self.a
        out = []
        if a.weighting in ("wfwhm", "nbstars"):
            out.append(f"{a.weighting} weighting (LSF weights differ from Siril's by design)")
        if a.rejection == "sigmedian" and a.weighting == "none":
            out.append("sigmedian rejection without weights")
        return out

    def load_pixels(self, path: str) -> np.ndarray:
        run_cmd(self.siril, "load", _q(path))
        return as_channels(self.siril.get_image_pixeldata())

    BREAKDOWN_FIELDS = ("n_diff", "max_abs_both_data", "n_one_zero", "n_diff_full", "max_abs_full",
                        "n_diff_partial", "max_abs_partial", "noise_lsf_adu", "noise_siril_adu", "noise_ratio",
                        "rows_noise_differ")

    def stage_parity(self, step) -> None:
        """LSF vs native Siril stacks, per channel. Same frames, same reference and settings, so the
        Stack of an lsf_init file is expected to equal Siril's stack exactly (see known_exceptions).
        Both images are read from the files in the same orientation; when they differ, a diff image
        <p>_diff_<pair>.fit (LSF - Siril) is written and the differences are broken down by coverage,
        with Siril's Noise Estimation computed on both images, row by row."""
        pairs = [(f"set_{k}", f"lsf_set_{k}.lsf.fits", f"siril_set_{k}.fit") for k in range(1, N_SETS + 1)]
        pairs += [(f"cum_{cum_tag(k)}", f"lsf_cum_{cum_tag(k)}.lsf.fits", f"siril_cum_{cum_tag(k)}.fit")
                  for k in range(1, N_SETS + 1)]
        exc = self.known_exceptions()
        if exc:
            self.log("Parity: known exceptions apply, differences are expected: " + "; ".join(exc), LogColor.SALMON)
        csv_path = self.p("parity.csv")
        mismatched = []
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["pair", "lsf_file", "siril_file", "channel", "max_abs_diff", "max_rel_diff"]
                       + list(self.BREAKDOWN_FIELDS) + ["known_exception"])
            for tag, lname, sname in pairs:
                lp, sp = self.p(lname), self.p(sname)
                if not (os.path.exists(lp) and os.path.exists(sp)):
                    self.log(f"  parity {tag}: missing {lname if not os.path.exists(lp) else sname}, skipped",
                             LogColor.SALMON)
                    continue
                step(f"Parity {tag}")
                try:
                    x, y = fits_image(lp), fits_image(sp)
                    if x is None or y is None:      # tile-compressed: Siril reads both the same way
                        x, y = self.load_pixels(lp), self.load_pixels(sp)
                        count = None
                    else:
                        count = fits_image(lp, "LSF_COUNT")
                    diffs = channel_diffs(x, y)
                    brk = parity_breakdown(x, y, count, header_int(lp, "STACKCNT"))
                    dpath = self.p(f"diff_{tag}.fit")
                    if any(b["n_diff"] for b in brk):
                        write_fits_float32(dpath, x - y, (f"LSF - Siril: {lname} - {sname}",))
                    elif os.path.exists(dpath):     # left over from an earlier run
                        os.remove(dpath)
                    del x, y, count
                except Exception as e:  # noqa: BLE001
                    self.log(f"  parity {tag}: failed ({e})", LogColor.SALMON)
                    w.writerow([tag, lname, sname, "", "", ""] + [""] * len(self.BREAKDOWN_FIELDS) + ["; ".join(exc)])
                    continue
                for c, ((ab, rel), b) in enumerate(zip(diffs, brk)):
                    b["noise_lsf_adu"], b["noise_siril_adu"] = b["noise_lsf"] * 65535, b["noise_siril"] * 65535
                    w.writerow([tag, lname, sname, c, repr(ab), repr(rel)]
                               + ["" if b.get(k) is None else repr(b[k]) for k in self.BREAKDOWN_FIELDS]
                               + ["; ".join(exc)])
                    self.log(f"  parity {tag:<9s} c{c}: max abs {ab:.3g} ({b['n_diff']} pixels differ"
                             + (f": {b['n_diff_full']} full / {b['n_diff_partial']} partial coverage"
                                if "n_diff_full" in b else "")
                             + f", {b['n_one_zero']} with data in only one)  Noise Estimation "
                             f"LSF {b['noise_lsf_adu']:.4f} / Siril {b['noise_siril_adu']:.4f} ADU "
                             f"(ratio {b['noise_ratio']:.6f})",
                             LogColor.GREEN if b["n_diff"] == 0 else LogColor.SALMON)
                    if b["worst_rows"]:
                        self.log("      rows whose noise differs most (file row: relative difference): "
                                 + ", ".join(f"{j + 1}: {r:.2e}" for j, r in b["worst_rows"]))
                    if b["n_diff"] and tag.startswith("cum_"):
                        mismatched.append(f"{tag} c{c}")
                fh.flush()
                gc.collect()
        self.log(f"Parity table written to {csv_path}", LogColor.GREEN)
        if mismatched:
            self.log("Full stacks that differ from Siril's: " + ", ".join(mismatched)
                     + (" (known exceptions apply)" if exc else " -- unexpected, see the diff_*.fit images"),
                     LogColor.SALMON if exc else LogColor.RED)

    SSIM_FILES = (("lsf_cum_1234", "lsf_cum_1234.lsf.fits"), ("lsf_4", "lsf_4.lsf.fits"),
                  ("lsf_sos", "lsf_sos.lsf.fits"), ("siril_cum_1234", "siril_cum_1234.fit"))
    SSIM_FIELDS = ("ssim", "one_minus_ssim", "ssim_p0_1", "ssim_min", "ms_ssim", "windows")

    def siril_image(self, path: str) -> tuple:
        """A tile-compressed file, read through Siril, in fits_memmap's form."""
        run_cmd(self.siril, "load", _q(path))
        a = np.asarray(self.siril.get_image_pixeldata(), dtype=np.float32)
        return (a[None] if a.ndim == 2 else a), 1.0, 0.0

    def stage_ssim(self, step) -> None:
        """SSIM and MS-SSIM, per channel, between every pair of the four full stacks (SSIM_FILES). The
        transform and data range come from siril_cum_1234 (else the first file found), so every pair
        is scored the same way. Files are read one channel at a time; when one is tile-compressed,
        all are read through Siril so that they share one orientation."""
        a = self.a
        files = [(lbl, saved_path(self.p(name))) for lbl, name in self.SSIM_FILES]
        for lbl, path in files:
            if not os.path.exists(path):
                self.log(f"  ssim: {os.path.basename(path)} missing, left out", LogColor.SALMON)
        files = [(lbl, path) for lbl, path in files if os.path.exists(path)]
        if len(files) < 2:
            self.log("ssim: fewer than two of the full stacks exist, skipped.", LogColor.SALMON)
            return
        via_siril = any(fits_memmap(path) is None for _, path in files)
        open_image = self.siril_image if via_siril else fits_memmap
        ref_label, ref_path = next(((lbl, p) for lbl, p in files if lbl == "siril_cum_1234"), files[0])
        step(f"SSIM reference transform from {ref_label}")
        ref = open_image(ref_path)
        transforms = []
        for c in range(ref[0].shape[0]):
            x = channel_f32(ref, c)
            transforms.append(ssim_reference(x[np.isfinite(x) & (x != 0)], a.ssim_domain))
        del ref, x
        self.log(f"SSIM in the {a.ssim_domain} domain, transform from {ref_label}: "
                 + "  ".join(f"c{c} range {L:.4g}" + ("" if c0 is None else f" shadows {c0:.4g} midtones {m:.4g}")
                             for c, (c0, m, L) in enumerate(transforms)))
        csv_path = self.p("ssim.csv")
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["pair", "file_a", "file_b", "channel", "domain", "data_range"] + list(self.SSIM_FIELDS)
                       + ["worst_blocks", "seconds"])
            for (la, pa), (lb, pb) in itertools.combinations(files, 2):
                tag = f"{la}_vs_{lb}"
                step(f"SSIM {tag}")
                names = [os.path.basename(pa), os.path.basename(pb)]
                try:
                    A, B = open_image(pa), open_image(pb)
                    if A[0].shape != B[0].shape:
                        raise RuntimeError(f"shape mismatch {A[0].shape} vs {B[0].shape}")
                    maps = []
                    for c in range(A[0].shape[0]):
                        t = time.perf_counter()
                        r = ssim_channel(channel_f32(A, c), channel_f32(B, c), transforms[min(c, len(transforms) - 1)],
                                         keep_map=a.ssim_maps)
                        r["one_minus_ssim"] = 1 - r["ssim"]
                        if a.ssim_maps:
                            maps.append(r.pop("dissim"))
                        dt = time.perf_counter() - t
                        worst = "; ".join(f"{y}:{x}={v:.6f}" for y, x, v in r["worst"])
                        w.writerow([tag] + names + [c, a.ssim_domain, repr(r["data_range"])]
                                   + [repr(r[k]) for k in self.SSIM_FIELDS] + [worst, f"{dt:.1f}"])
                        self.log(f"  ssim {tag} c{c}: SSIM {r['ssim']:.6f} (1-SSIM {r['one_minus_ssim']:.3g}, "
                                 f"0.1% {r['ssim_p0_1']:.4f}, min {r['ssim_min']:.4f})  MS-SSIM {r['ms_ssim']:.6f}"
                                 f"  [{dt:.1f} s]")
                        if r["worst"] and r["ssim_min"] < 1:
                            self.log("      lowest 64x64 blocks (file row:column = mean SSIM): " + worst)
                    del A, B
                    if maps:
                        write_fits_float32(self.p(f"ssim_{tag}.fit"), np.stack(maps),
                                           (f"1 - SSIM ({a.ssim_domain}): {names[0]} vs {names[1]}",))
                        del maps
                except Exception as e:  # noqa: BLE001
                    self.log(f"  ssim {tag}: failed ({e})", LogColor.SALMON)
                    w.writerow([tag] + names + ["", a.ssim_domain, ""] + [""] * len(self.SSIM_FIELDS) + ["", ""])
                fh.flush()
                gc.collect()
        self.log(f"SSIM table written to {csv_path}", LogColor.GREEN)

    # ---- behaviour tests --------------------------------------------------------
    def try_cmd(self, *parts) -> bool:
        """A command that may be refused: True if Siril reports success. The tests below judge the
        outcome from the files, so a refusal reported either way is fine."""
        try:
            run_cmd(self.siril, *parts)
            return True
        except RuntimeError:
            return False

    def stage_flip(self, step) -> None:
        """An LSF file whose Stack is flipped by an outside tool that keeps the old
        DATASUM. lsf_append must refuse it, lsf_util -repair must restore the Stack bit-exactly, and an
        edit of one pixel kept with -adopt must leave the rest of the Stack unchanged."""
        src = next((self.p(n) for n in ("lsf_cum_1234.lsf.fits", "lsf_set_1.lsf.fits", "lsf_4.lsf.fits")
                    if os.path.exists(self.p(n))), None)
        if src is None:
            self.log("flip: no LSF product to test with, skipped.", LogColor.SALMON)
            return
        step("Flip, repair and adopt test")
        rows, failed = [], []

        def chk(name, ok, detail=""):
            rows.append(("flip", name, "ok" if ok else "FAIL", detail))
            if not ok:
                failed.append(f"{name} ({detail})")

        work = os.path.join(self.workdir, "flip.lsf.fits")
        shutil.copy2(src, work)
        orig = primary_bytes(work)
        gen0 = header_int(work, "LSFGEN")
        edit_primary(work, lambda a: a[:, ::-1, :])
        chk("flipped_differs", primary_bytes(work) != orig, os.path.basename(src))
        self.try_cmd("lsf_util", _q(work), "-check")
        self.try_cmd("lsf_append", _q(work), _q(self.seqname), "-noreg")
        chk("append_refused", header_int(work, "LSFGEN") == gen0, f"LSFGEN {header_int(work, 'LSFGEN')}")
        self.try_cmd("lsf_util", _q(work), "-repair")
        chk("repair_bit_exact", primary_bytes(work) == orig, "primary data compared byte by byte")
        chk("repair_new_generation", header_int(work, "LSFGEN") == gen0 + 1, f"LSFGEN {header_int(work, 'LSFGEN')}")
        chk("passes_after_repair", all(ok for _, ok, _ in check_lsf_file(work)[0]), "structure and LSFHASH")

        # adopt: one pixel with data, set to a value strictly inside the range of the Stack
        x = fits_image(work)
        ch, ys, xs = np.nonzero((x > np.quantile(x, 0.4)) & (x < np.quantile(x, 0.6)))
        c, y, xx = int(ch[len(ch) // 2]), int(ys[len(ys) // 2]), int(xs[len(xs) // 2])
        new_value = np.float32((x[c, y, xx] + np.quantile(x, 0.6)) / 2)

        def poke(a):
            a[c, y, xx] = new_value
            return a
        edit_primary(work, poke)
        before = fits_image(work)
        self.try_cmd("lsf_util", _q(work), "-adopt")
        after = fits_image(work)
        changed = np.argwhere(after != before)
        chk("adopt_keeps_the_edit", after[c, y, xx] == new_value, f"pixel ({c},{y},{xx})")
        chk("adopt_changes_nothing_else", len(changed) == 0, f"{len(changed)} other pixels differ")
        chk("adopt_new_generation", header_int(work, "LSFGEN") == gen0 + 2, f"LSFGEN {header_int(work, 'LSFGEN')}")
        chk("passes_after_adopt", all(ok for _, ok, _ in check_lsf_file(work)[0]), "structure and LSFHASH")
        self.write_test_csv("flip", rows, failed)

    def write_test_csv(self, name: str, rows: list, failed: list) -> None:
        csv_path = self.p(f"{name}.csv")
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["test", "check", "result", "detail"])
            w.writerows(rows)
        if failed:
            self.log(f"{name} test FAILED: " + "; ".join(failed) + f" -> {csv_path}", LogColor.RED)
            raise RuntimeError(f"{name} test failed: " + "; ".join(failed))
        self.log(f"{name} test passed -> {csv_path}", LogColor.GREEN)

    # ---- noise ----------------------------------------------------------------
    def measure_own_noise(self, channels: int, shape) -> list:
        """Independent noise estimate from a few bands of rows (never the whole image)."""
        if not self._own_noise_ok:
            return [None] * channels
        try:
            _, h, w = shape
            band_h = max(2, min(64, h))
            ys = sorted({int(y) for y in np.linspace(0, max(h - band_h, 0), 8)})
            bands = [np.asarray(self.siril.get_image_pixeldata(shape=[0, y, w, band_h])) for y in ys]
            return own_noise_from_bands(bands, channels)
        except Exception as e:  # noqa: BLE001
            self._own_noise_ok = False
            self.log(f"Independent noise estimate unavailable ({e}); continuing with Siril's statistics only.",
                     LogColor.SALMON)
            return [None] * channels

    def measure_file(self, label, kind, path, expected, writer, fh) -> list:
        """Load one product (LSF files directly), run `stat`, append rows to the CSV."""
        a = self.a
        stackcnt = header_int(path, "STACKCNT")
        run_cmd(self.siril, "load", _q(path))
        shape = self.siril.get_image_shape()
        channels = int(shape[0]) if shape else 1
        noise_est = [None] * channels
        if a.cfa or a.stat_main:
            # `stat -cfa` works on the CFA sub-channels and `stat main` has no noise estimate:
            # a plain `stat` first gives the value Noise Estimation reports
            run_cmd(self.siril, "stat")
            noise_est = self.noise_estimation_adu(channels, path)
        run_cmd(self.siril, "stat", *(["-cfa"] if a.cfa else []), *(["main"] if a.stat_main else []))
        if not (a.cfa or a.stat_main):
            noise_est = self.noise_estimation_adu(channels, path)
        known = [v for v in noise_est if v is not None]
        noise_mean = sum(known) / len(known) if len(known) == channels else None
        own = self.measure_own_noise(channels, shape) if (a.own_noise and shape) else [None] * channels
        summary = []
        for c in range(channels):
            st = None
            try:
                st = self.siril.get_image_stats(c)
            except Exception as e:  # noqa: BLE001
                self.log(f"  {label}: statistics for channel {c} unavailable ({e})", LogColor.SALMON)
            vals = {k: (getattr(st, k, None) if st is not None else None) for k in STAT_FIELDS}
            med, on = vals["median"], own[c]
            writer.writerow([label, kind, os.path.basename(path), expected if expected is not None else "",
                             stackcnt if stackcnt is not None else "", c]
                            + ["" if vals[k] is None or float(vals[k]) == -999999.0 else repr(float(vals[k])) for k in STAT_FIELDS]
                            + ["" if noise_est[c] is None else repr(noise_est[c]),
                               "" if noise_mean is None else repr(noise_mean)]
                            + ["" if on is None else repr(on),
                               "" if (on is None or not med) else repr(on / float(med))])
            summary.append((noise_est[c], on))
        fh.flush()
        return summary

    def noise_estimation_adu(self, channels: int, path: str) -> list:
        """Per-channel value of Tools > Image Analysis > Noise Estimation (the bgnoise command), in
        ADU as Siril logs it: x65535 for float images.

        Siril's statistics of the loaded image come first. In this build, after `stat` they come
        back from the Python API as -999999 (not computed), so the fallback is
        analyse_image_from_file(), which runs the same estimator (STATS_SIGMEAN, FnNoise1) but only
        on channel 0, or green for colour. Any other channel uses noise_estimation(), a port of
        FnNoise1 that matches Siril's value to 4 decimals, on the pixels of the file."""
        try:
            img = self.siril.get_image(with_pixels=False)
            is_float = int(getattr(img, "bitpix", -32)) in (-32, -64)
        except Exception:  # noqa: BLE001
            is_float = True
        scale = 65535.0 if is_float else 1.0
        out = []
        for c in range(channels):
            try:
                bg = getattr(self.siril.get_image_stats(c), "bgnoise", None)
            except Exception:  # noqa: BLE001
                bg = None
            out.append(float(bg) * scale if bg is not None and bg >= 0 else None)
        if all(v is not None for v in out):
            return out
        try:        # analyse_image_from_file returns bgnoise in [0, 1] for every image type
            bg = float(self.siril.analyse_image_from_file(path).bgnoise)
            ch = 0 if channels == 1 else 1
            if out[ch] is None and bg > 0:
                out[ch] = bg * 65535.0
        except Exception as e:  # noqa: BLE001
            self.log(f"  analyse_image_from_file failed ({e}); using the Python noise estimate.", LogColor.SALMON)
        if any(v is None for v in out):
            x = fits_image(path)
            for c in range(channels):
                if out[c] is None and x is not None and c < x.shape[0]:
                    out[c] = noise_estimation(x[c])[0] * 65535.0
        return out

    def stage_stats(self, step) -> None:
        csv_path = self.p("noise.csv")
        rows_seen = []
        files = [pr for pr in self.products if os.path.exists(pr[2])]
        have = {pr[2] for pr in files}
        for label, kind, name in self.known_outputs():      # products of skipped stages already on disk
            path = self.p(name)
            if path not in have and os.path.exists(path):
                files.append((label, kind, path, None))
        order = {name: i for i, (_, _, name) in enumerate(self.known_outputs())}
        files.sort(key=lambda pr: order.get(os.path.basename(pr[2])[len(self.a.prefix) + 1:], 99))
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["label", "kind", "file", "frames_expected", "stackcnt_header", "channel"]
                            + list(STAT_FIELDS) + ["noise_est_adu", "noise_est_mean_adu", "own_noise",
                                                     "own_noise_over_median"])
            for label, kind, path, expected in files:
                step(f"Noise statistics: {label}")
                try:
                    rows_seen.append((label, self.measure_file(label, kind, path, expected, writer, fh)))
                except Exception as e:  # noqa: BLE001
                    self.log(f"  {label}: could not be measured ({e})", LogColor.SALMON)
                    writer.writerow([label, kind, os.path.basename(path), "", "", "", *[""] * len(STAT_FIELDS),
                                     "", "", "", ""])
                gc.collect()
        self.log(f"Noise table written to {csv_path}", LogColor.GREEN)
        self.log("Noise summary, ADU (Siril Noise Estimation | independent estimate x65535, per channel):")
        fmt = lambda v, k=1.0: "n/a" if v is None else f"{float(v) * k:.6g}"        # noqa: E731
        for label, summ in rows_seen:
            self.log(f"  {label:<14s} " + "  ".join(f"[{fmt(b)} | {fmt(o, 65535.0)}]" for b, o in summ))

    @staticmethod
    def known_outputs() -> list:
        """(label, kind, file name without prefix) of every possible product, in report order."""
        tags = [cum_tag(k) for k in range(1, N_SETS + 1)]
        out = [(f"lsf_cum_{t}", "lsf_base", f"lsf_cum_{t}.lsf.fits") for t in tags]
        out += [(f"lsf_set_{k}", "lsf_base", f"lsf_set_{k}.lsf.fits") for k in range(1, N_SETS + 1)]
        out.append(("lsf_sos", "lsf_base", "lsf_sos.lsf.fits"))
        out += [(f"lsf_{k}", "lsf_iter", f"lsf_{k}.lsf.fits") for k in range(1, N_SETS + 1)]
        out += [(f"siril_cum_{t}", "siril", f"siril_cum_{t}.fit") for t in tags]
        out += [(f"siril_set_{k}", "siril", f"siril_set_{k}.fit") for k in range(1, N_SETS + 1)]
        return out

    # ---- driver -------------------------------------------------------------------
    def prepare(self) -> None:
        a, siril = self.a, self.siril
        if a.seq:
            run_cmd(siril, "load_seq", _q(a.seq))
        seq = self.current_seq()
        if seq is None:
            raise RuntimeError("No sequence is loaded in Siril (pass --seq or load one first).")
        self.seqname, self.nframes = seq.seqname, int(seq.number)
        if a.rgb_equal and int(getattr(seq, "nb_layers", 1)) != 3:
            a.rgb_equal = False
            self.log("Sequence is not 3-channel: RGB equalization turned off.")
        self.original_selected = [i for i in range(self.nframes) if seq.imgparam[i].incl]
        self.sets = split_sets(self.original_selected, N_SETS, a.split)
        self.wd = siril.get_siril_wd()
        self.outdir = os.path.abspath(a.outdir or os.path.join(self.wd, "lsf_compare"))

    def prime_stats_cache(self) -> None:
        """Make every stack of this run normalize with the same frame statistics. Siril keeps the
        statistics it computes in the .seq file, with 6 significant digits: the first command that
        needs them on a fresh sequence uses the full-precision values, later ones read the rounded
        cache. lsf_init and stack would then differ by float rounding and by the odd rejection
        decision. `seqstat ... full` computes and caches the normalization statistics (location and
        scale, IKSS and lite) of the selected frames, and reloading the sequence puts the cached
        values in memory for `stack` too; lsf_init reads them from the .seq file."""
        t = time.perf_counter()
        self.ensure_seq()
        run_cmd(self.siril, "seqstat", _q(self.seqname), _q(os.path.join(self.workdir, "seqstat.csv")), "full")
        run_cmd(self.siril, "load_seq", _q(self.seqname))
        self.log(f"Frame statistics cached in {self.seqname}.seq ({time.perf_counter() - t:.1f} s), so that "
                 f"LSF and Siril stacks use the same values.")

    def describe_plan(self) -> None:
        a = self.a
        self.log(f"Sequence '{self.seqname}': {len(self.original_selected)} of {self.nframes} frames selected, "
                 f"split {a.split} into sets of " + ", ".join(str(len(x)) for x in self.sets) + ".")
        for k, grp in enumerate(self.sets, 1):
            shown = ", ".join(f"{x + 1}-{y + 1}" if x != y else f"{x + 1}" for x, y in runs_of(grp))
            self.log(f"  set {k} ({len(grp)}): frames {shown if len(shown) < 200 else shown[:200] + ' ...'}")
        self.log(f"Output folder: {self.outdir}")
        self.log("Stages: " + ", ".join(st for st in STAGES if st not in a.skip))
        self.log("  lsf_init:   " + " ".join(self.lsf_init_cmd(self.p("lsf_set_1.lsf.fits"))))
        self.log("  lsf_append: " + " ".join(self.lsf_append_cmd(self.p("lsf_2.lsf.fits"))))
        self.log("  stack:      " + " ".join(self.siril_stack_cmd(self.p("siril_set_1.fit"))))
        for e in self.known_exceptions():
            self.log(f"  parity known exception: {e}", LogColor.SALMON)
        if "sos" not in a.skip and "base" in a.skip:
            self.log("  sos needs existing lsf_set_k files from an earlier run.")
        if "parity" not in a.skip and ("base" in a.skip or "siril" in a.skip):
            self.log("  parity uses existing files for the skipped base/siril stage.")
        if min(len(x) for x in self.sets) < 5:
            self.log("WARNING: sets of fewer than 5 frames: sigma-type rejection does almost nothing there.",
                     LogColor.SALMON)
        existing = [saved_path(f) for f in self.planned_files() if os.path.exists(saved_path(f))]
        if existing and not a.overwrite:
            raise RuntimeError("These products already exist (use --overwrite to replace them): "
                               + ", ".join(os.path.basename(f) for f in existing))

    def run(self) -> None:
        a, siril = self.a, self.siril
        ok = False
        try:
            self.prepare()
            self.describe_plan()
            if a.dry_run:
                self.log("Dry run: planned files: " + ", ".join(os.path.basename(f) for f in self.planned_files()))
                ok = True
                return
            os.makedirs(self.outdir, exist_ok=True)
            self.workdir = tempfile.mkdtemp(prefix=".lsfcmp_", dir=self.outdir)
            if "base" not in a.skip or "siril" not in a.skip:
                self.prime_stats_cache()

            weights = {"base": 7, "siril": 7, "sos": 1, "iter": 4, "check": 1, "flip": 1, "parity": 8, "ssim": 7,
                       "stats": 21}
            stages = {"base": self.stage_base, "siril": self.stage_siril, "sos": self.stage_sos,
                      "iter": self.stage_iter, "check": self.stage_check, "parity": self.stage_parity,
                      "ssim": self.stage_ssim, "stats": self.stage_stats, "flip": self.stage_flip}
            todo = [st for st in STAGES if st not in a.skip]
            total_steps = sum(weights[st] for st in todo)
            done = [0]

            def step(msg):
                siril.update_progress(msg, min(0.99, done[0] / max(total_steps, 1)))
                done[0] += 1

            for st in todo:
                stages[st](step)
            siril.update_progress("Done", 1.0)
            self.log(f"All done in {time.perf_counter() - self.t0:.0f} s. Products are in {self.outdir}",
                     LogColor.GREEN)
            ok = True
        finally:
            self.cleanup(ok)

    def cleanup(self, ok: bool) -> None:
        """Always: scratch files gone, original selection back, memory released."""
        try:
            if self.workdir:
                shutil.rmtree(self.workdir, ignore_errors=True)
            if self.selection_touched and self.seqname and self.original_selected:
                try:
                    self.apply_selection(self.original_selected)      # reloads the sequence if needed
                except Exception as e:  # noqa: BLE001
                    self.log(f"Could not restore your original frame selection automatically: {e}. "
                             f"It was frames {runs_of(self.original_selected)} (0-based).", LogColor.SALMON)
            if not ok:
                self.log(f"Run stopped early; scratch files were removed. Finished products remain in {self.outdir}.",
                         LogColor.SALMON)
        finally:
            gc.collect()


def confirm_gui(cmp_: Comparison) -> bool:
    """Scripts-menu launch with no arguments: show the plan once and ask before a long run."""
    try:
        import tkinter as tk
        from tkinter import messagebox
        cmp_.prepare()
        root = tk.Tk()
        root.withdraw()
        msg = (f"Sequence '{cmp_.seqname}': {len(cmp_.original_selected)} selected frames, sets of "
               + ", ".join(str(len(x)) for x in cmp_.sets)
               + f".\n\nThis will build LSF and Siril per-set and cumulative stacks, a stack of stacks and 4 "
                 f"iterative LSF files, check them and export parity, SSIM and noise tables to:\n{cmp_.outdir}\n\n"
                 f"This can take a long time. Continue?")
        ok = messagebox.askokcancel("LSF compare", msg)
        root.destroy()
        return bool(ok)
    except Exception:  # noqa: BLE001
        return True


def main():
    siril = s.SirilInterface()
    siril.connect()
    try:
        argv = sys.argv[1:]
        args = build_arg_parser().parse_args(argv)
        job = Comparison(siril, args)
        if not argv and not siril.is_cli() and not confirm_gui(job):
            siril.log("LSF compare cancelled by user.")
            return
        job.run()
    except Exception as e:  # noqa: BLE001
        siril.log(f"LSF compare failed: {e}", LogColor.RED)
        traceback.print_exc()
        try:
            siril.error_messagebox(f"LSF compare failed:\n{e}")
        except Exception:  # noqa: BLE001
            pass
        sys.exit(1)
    finally:
        siril.disconnect()


if __name__ == "__main__":
    main()
