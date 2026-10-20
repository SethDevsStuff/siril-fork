#!/usr/bin/env python3
"""
lsf_compare.py -- end-to-end test of Siril's native Living Stack Format (LSF v0.4) commands
(lsf_init, lsf_append, lsf_export): fresh LSF stacks versus iterative appends, parity against
Siril's own `stack` command, v0.4 structural checks, and a noise estimate for every product.

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
                                    LSF_ACC NAXIS3 == 4*P, LSFUID/LSFGEN identical across LSF HDUs,
                                    LSFGEN +1 along lsf_1..lsf_4, STACKCNT == ledger rows,
                                    LSFHASH == sha256(primary data), plus `fitsverify` if on PATH.
                                    Raises if any check fails.
Stage parity <p>_parity.csv         lsf_set_k vs siril_set_k and lsf_cum_X vs siril_cum_X:
                                    per-channel max abs and max rel difference. Known, expected
                                    exceptions (logged): wfwhm / nbstars weighting, and sigmedian
                                    rejection without weights.
Stage stats  <p>_noise.csv          one row per file and channel; the .lsf.fits files are loaded
                                    directly (this also tests that Siril opens one as an image).

--export-fits additionally writes <name>.fit next to every LSF product with `lsf_export`.

The stack of set 1 is the same frames with the same settings whichever way it is stacked, so it
is computed once (lsf_set_1) and copied to lsf_cum_1 and lsf_1 (siril_set_1 -> siril_cum_1).

NOISE
-----
For every product the script loads it in Siril, runs `stat [-cfa] [main]` and reads the per-channel
statistics back through the Python API. Siril's background noise estimate is the `bgnoise` column.
The script also measures an independent per-pixel noise (robust MAD of neighbour differences on a
few horizontal bands; --no-own-noise to skip).

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
import math
import os
import re
import shutil
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
STAGES = ("base", "siril", "sos", "iter", "check", "parity", "stats")   # run order
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
    chk("acc_naxis3_4P", n3 == 4 * p, f"NAXIS3={n3} P={p}")
    planes = [acc.get(f"LSFPL{i}") for i in range(1, n3 + 1)]
    chk("acc_planes_named", all(planes), ",".join(str(x) for x in planes))
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


def run_fitsverify(path: str):
    """(ok, detail) from fitsverify, or None if it is not on PATH."""
    exe = shutil.which("fitsverify")
    if not exe:
        return None
    r = subprocess.run([exe, "-q", path], capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip()
    m = re.search(r"(\d+)\s+error", out)
    errors = int(m.group(1)) if m else (0 if r.returncode == 0 else 1)
    return errors == 0, out.splitlines()[-1] if out else f"exit {r.returncode}"


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
                        "parity (LSF vs Siril), stats (noise export).")
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
        for st, name in (("check", "check.csv"), ("parity", "parity.csv"), ("stats", "noise.csv")):
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
        if not os.path.exists(path):
            raise RuntimeError(f"{label}: the command reported success but {os.path.basename(path)} is missing.")
        self.check_count(path, expected, label)
        self.products.append((label, kind, path, expected))
        self.done[label] = path
        if kind.startswith("lsf") and self.a.export_fits:
            self.export_plain(path)

    def export_plain(self, lsf_path: str) -> None:
        """lsf_export of an LSF file to <name>.fit. A failed export is a warning."""
        dest = plain_path(lsf_path)
        try:
            if os.path.exists(dest):
                os.remove(dest)
            run_cmd(self.siril, "lsf_export", _q(lsf_path), _q(dest))
            self.log(f"  lsf_export -> {os.path.basename(dest)}")
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
            if os.path.exists(dest):
                os.remove(dest)
            self.apply_selection(S[k - 1])
            self.timed(self.siril_stack_cmd(dest), f"siril_set_{k}", f"{len(S[k - 1])} frames stacked")
            self.register(f"siril_set_{k}", "siril", dest, len(S[k - 1]))
        self.copy_product("siril_set_1", self.p("siril_cum_1.fit"), "siril_cum_1", "siril", len(S[0]),
                          "the running stack over set 1 IS the stack of set 1")
        for k in range(2, N_SETS + 1):
            step(f"Siril stack of sets 1-{k}")
            frames = sorted(i for grp in S[:k] for i in grp)
            dest = self.p(f"siril_cum_{cum_tag(k)}.fit")
            if os.path.exists(dest):
                os.remove(dest)
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

    def stage_parity(self, step) -> None:
        """LSF vs native Siril stacks, per channel."""
        pairs = [(f"set_{k}", f"lsf_set_{k}.lsf.fits", f"siril_set_{k}.fit") for k in range(1, N_SETS + 1)]
        pairs += [(f"cum_{cum_tag(k)}", f"lsf_cum_{cum_tag(k)}.lsf.fits", f"siril_cum_{cum_tag(k)}.fit")
                  for k in range(1, N_SETS + 1)]
        exc = self.known_exceptions()
        if exc:
            self.log("Parity: known exceptions apply, differences are expected: " + "; ".join(exc), LogColor.SALMON)
        csv_path = self.p("parity.csv")
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["pair", "lsf_file", "siril_file", "channel", "max_abs_diff", "max_rel_diff",
                        "known_exception"])
            for tag, lname, sname in pairs:
                lp, sp = self.p(lname), self.p(sname)
                if not (os.path.exists(lp) and os.path.exists(sp)):
                    self.log(f"  parity {tag}: missing {lname if not os.path.exists(lp) else sname}, skipped",
                             LogColor.SALMON)
                    continue
                step(f"Parity {tag}")
                try:
                    x = self.load_pixels(lp)
                    y = self.load_pixels(sp)
                    diffs = channel_diffs(x, y)
                    del x, y
                except Exception as e:  # noqa: BLE001
                    self.log(f"  parity {tag}: failed ({e})", LogColor.SALMON)
                    w.writerow([tag, lname, sname, "", "", "", "; ".join(exc)])
                    continue
                for c, (ab, rel) in enumerate(diffs):
                    w.writerow([tag, lname, sname, c, repr(ab), repr(rel), "; ".join(exc)])
                self.log(f"  parity {tag:<9s} " + "  ".join(f"[c{c}: abs {ab:.3g} rel {rel:.3g}]"
                                                         for c, (ab, rel) in enumerate(diffs)))
                fh.flush()
                gc.collect()
        self.log(f"Parity table written to {csv_path}", LogColor.GREEN)

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
        run_cmd(self.siril, "stat", *(["-cfa"] if a.cfa else []), *(["main"] if a.stat_main else []))
        shape = self.siril.get_image_shape()
        channels = int(shape[0]) if shape else 1
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
                            + ["" if vals[k] is None else repr(float(vals[k])) for k in STAT_FIELDS]
                            + ["" if on is None else repr(on),
                               "" if (on is None or not med) else repr(on / float(med))])
            summary.append((vals["bgnoise"], on))
        fh.flush()
        return summary

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
                            + list(STAT_FIELDS) + ["own_noise", "own_noise_over_median"])
            for label, kind, path, expected in files:
                step(f"Noise statistics: {label}")
                try:
                    rows_seen.append((label, self.measure_file(label, kind, path, expected, writer, fh)))
                except Exception as e:  # noqa: BLE001
                    self.log(f"  {label}: could not be measured ({e})", LogColor.SALMON)
                    writer.writerow([label, kind, os.path.basename(path), "", "", "", *[""] * len(STAT_FIELDS),
                                     "", ""])
                gc.collect()
        self.log(f"Noise table written to {csv_path}", LogColor.GREEN)
        self.log("Noise summary (Siril bgnoise | independent estimate, per channel):")
        fmt = lambda v: "n/a" if v is None else f"{float(v):.6g}"        # noqa: E731
        for label, summ in rows_seen:
            self.log(f"  {label:<14s} " + "  ".join(f"[{fmt(b)} | {fmt(o)}]" for b, o in summ))

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
        existing = [f for f in self.planned_files() if os.path.exists(f)]
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

            weights = {"base": 7, "siril": 7, "sos": 1, "iter": 4, "check": 1, "parity": 8, "stats": 21}
            stages = {"base": self.stage_base, "siril": self.stage_siril, "sos": self.stage_sos,
                      "iter": self.stage_iter, "check": self.stage_check, "parity": self.stage_parity,
                      "stats": self.stage_stats}
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
                 f"iterative LSF files, check them and export parity and noise tables to:\n{cmp_.outdir}\n\n"
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
