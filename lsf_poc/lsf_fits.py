"""
lsf_fits.py -- Living Stack Format (LSF) FITS profile (live_stack_fits.md, v0.3).

Container layer only. The LSF data model (accumulator, frame ledger, canonical grid,
normalization, settings) is unchanged; this module maps it onto the HDU layout of the
FITS profile and back:

    HDU 0  PRIMARY      Stack (weighted mean, float32, uncompressed, WCS)
    HDU 1  LSF_WEIGHT   sum of weights, float32
    HDU 2  LSF_M2       Welford M2, float32
    HDU 3  LSF_W2       sum of squared weights, float32 (always per channel)
    HDU 4  LSF_COUNT    contribution count, uint16 (BZERO 32768), saturating
    HDU 5  LSF_SETTINGS BINTABLE  KEY / TYPE / VALUE
           ...          unrecognized HDUs from the file being rewritten (C7)
    HDU n  LSF_LEDGER   BINTABLE, always last, one row per contributing frame

The rest of the code base works on an in-memory model (not on HDUs): read_lsf_fits() returns
{"images": channel-last arrays, "properties": the JSON blocks LSF:FrameLedger / CanonicalGrid /
Settings / OutputNorm / NormBaseline, "meta": what a rewrite must carry forward (LSFUID,
LSFGEN, unknown keywords, unknown HDUs, unknown ledger columns, ...)}, and write_lsf_fits()
takes the same images and properties plus that meta.

Rules from the spec that this module enforces:
  C1  HDU 0 is the real, uncompressed float32 Stack.
  C2  Only standard IMAGE / BINTABLE HDUs, no HIERARCH needed.
  C5  LSF keywords <= 8 characters; extension names use the LSF_ prefix.
  C6  HDUs are located by EXTNAME, never by position.
  C7  Unknown HDUs, primary keywords, settings rows and ledger columns survive a rewrite.
  S4  No BAYERPAT / XBAYROFF / YBAYROFF in the primary; ROWORDER in every image HDU;
      LSFNMIN / LSFNMAX written so they round-trip exactly.
  S8  Pre-merge validation (checks 1-6) via validate_hdul() / validate_lsf_fits().

Requires astropy.
"""

from __future__ import annotations

import datetime
import json
import math
import os
import re
import uuid

import numpy as np

try:
    from astropy.io import fits
    ASTROPY_IMPORT_ERROR = None
except ImportError as _e:       # reported by _require_astropy() when a FITS function is used
    fits = None
    ASTROPY_IMPORT_ERROR = _e

PROFILE_VERSION = "0.3"
SUPPORTED_MAJOR = 0
DEFAULT_ROW_ORDER = "TOP-DOWN"

EXT_WEIGHT, EXT_M2, EXT_W2, EXT_COUNT = "LSF_WEIGHT", "LSF_M2", "LSF_W2", "LSF_COUNT"
EXT_SETTINGS, EXT_LEDGER = "LSF_SETTINGS", "LSF_LEDGER"
# in-memory image id -> EXTNAME
# (dict order is the HDU order of Section 3: WEIGHT, M2, W2, COUNT)
IMAGE_EXTS = {"WeightSum": EXT_WEIGHT, "SumSq": EXT_M2, "SumW2": EXT_W2, "ContributionCount": EXT_COUNT}
KNOWN_EXTS = (EXT_WEIGHT, EXT_M2, EXT_W2, EXT_COUNT, EXT_SETTINGS, EXT_LEDGER)

COMPRESS_MODES = ("gzip2",)       # float planes GZIP_2 (lossless), count plane RICE_1

INT_NULL = -2147483648            # TNULL for the J columns of the ledger

# Keywords owned by the writer: never copied from a caller or from the old file.
_STRUCTURAL_RE = re.compile(
    r"^(SIMPLE|BITPIX|NAXIS\d*|EXTEND|XTENSION|PCOUNT|GCOUNT|BSCALE|BZERO|BLANK|CHECKSUM|DATASUM"
    r"|EXTNAME|EXTVER|ROWORDER|TFIELDS|THEAP|T(FORM|TYPE|UNIT|NULL|SCAL|ZERO|DISP|DIM)\d+|Z[A-Z]+\d*)$")
LSF_ID_KEYS = ("LSF", "LSFVER", "LSFUID", "LSFGEN", "LSFSENS", "LSFBAYR", "LSFDMOD",
               "LSFNORM", "LSFNMIN", "LSFNMAX", "LSFDSUM", "LSFBUSY")
BAYER_KEYS = ("BAYERPAT", "XBAYROFF", "YBAYROFF")
_COMMENTARY = ("HISTORY", "COMMENT", "")
_WCS_RE = re.compile(
    r"^((CTYPE|CRVAL|CRPIX|CRDER|CSYER|CUNIT|CDELT|CROTA)[1-3]|CD[1-3]_[1-3]|PC[1-3]_[1-3]"
    r"|EQUINOX|EPOCH|RADESYS|LONPOLE|LATPOLE|WCSAXES|(A|B|AP|BP)_(\d+_\d+|ORDER|DMAX))$")

_BAYER_PATTERNS = ("RGGB", "GRBG", "GBRG", "BGGR")
_DMOD_FROM_METHOD = {"mono": "NONE", "debayer_first": "DEBAYER", "cfa_drizzle": "CFADRIZ"}
_METHOD_FROM_DMOD = {v: k for k, v in _DMOD_FROM_METHOD.items()}


class LsfFitsError(RuntimeError):
    pass


class LsfValidationError(LsfFitsError):
    """A Section 8 pre-merge check failed. `check` is the check number (1-6)."""

    def __init__(self, check: int, message: str):
        self.check, self.detail = check, message
        super().__init__(f"LSF validation failed (check {check}): {message}. "
                         f"The file was left untouched; it is still readable as an image.")


def _require_astropy():
    if fits is None:
        raise LsfFitsError("The FITS container needs astropy (pip install astropy). "
                           f"Import error: {ASTROPY_IMPORT_ERROR}")


def is_fits_file(path: str) -> bool:
    try:
        with open(path, "rb") as f:
            return f.read(9) == b"SIMPLE  ="
    except OSError:
        return False


def is_wcs_key(key: str) -> bool:
    return bool(_WCS_RE.match(key.upper()))


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def _now_utc() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _py(v):
    """numpy scalar -> python scalar (astropy rejects some numpy types in headers)."""
    return v.item() if isinstance(v, np.generic) else v


def _exact_float_card(key: str, value: float, comment: str = ""):
    """A real-valued card that round-trips exactly. astropy's own formatter keeps 16
    significant digits, one short of what an arbitrary float64 needs; repr() is the
    shortest string that parses back to the same double."""
    s = repr(float(value)).upper()
    if "." not in s and "E" not in s:
        s += ".0"
    text = f"{key:<8}= {s:>20}" + (f" / {comment}" if comment else "")
    return fits.Card.fromstring(text[:80])


def _set_card(hdr, key, value, comment=""):
    value = _py(value)
    if isinstance(value, float) and math.isfinite(value) and key in ("LSFNMIN", "LSFNMAX"):
        hdr.append(_exact_float_card(key, value, comment), end=True)
    else:
        hdr[key] = (value, comment)


def datasum_of_planes(arr: np.ndarray, chunk_elems: int = 1 << 22) -> int:
    """FITS DATASUM (32-bit ones-complement sum of the big-endian 32-bit words of the data
    unit) of a float32 / int32 image, computed in chunks so a large plane is never copied
    whole. Word order does not matter for the sum, so a channel-last array gives the same
    value as the channel-first array that is actually written."""
    if arr.dtype.itemsize != 4:
        raise LsfFitsError("datasum_of_planes only handles 4-byte samples")
    total = 0
    flat_rows = arr.reshape(-1, arr.shape[-1]) if arr.ndim > 1 else arr.reshape(1, -1)
    rows_per = max(1, chunk_elems // max(1, flat_rows.shape[1]))
    for r0 in range(0, flat_rows.shape[0], rows_per):
        blk = np.ascontiguousarray(flat_rows[r0:r0 + rows_per], dtype=arr.dtype.newbyteorder(">"))
        total += int(blk.view(">u4").sum(dtype=np.uint64))
    s = total % 0xFFFFFFFF
    return 0xFFFFFFFF if (s == 0 and total != 0) else s


def _to_fits_planes(arr: np.ndarray) -> np.ndarray:
    """In-memory channel-last (H, W[, C]) -> FITS channel-first contiguous (C, H, W) / (H, W)."""
    if arr.ndim == 3:
        return np.ascontiguousarray(np.transpose(arr, (2, 0, 1)))
    return np.ascontiguousarray(arr)


def _from_fits_planes(data) -> np.ndarray:
    """FITS (C, H, W) / (H, W) -> channel-last view (no copy)."""
    return np.transpose(data, (1, 2, 0)) if data.ndim == 3 else data


def native(arr) -> np.ndarray:
    """C-contiguous native-byte-order copy of an array that may be a big-endian FITS memmap."""
    a = np.asarray(arr)
    return np.ascontiguousarray(a, dtype=a.dtype.newbyteorder("="))


def _hdu_name(h) -> str:
    return str(h.header.get("EXTNAME", "") or "").strip().upper()


def _find_hdus(hdul) -> dict:
    """EXTNAME -> HDU, for every extension (C6: never by position)."""
    out = {}
    for i, h in enumerate(hdul):
        if i == 0:
            continue
        out.setdefault(_hdu_name(h), h)
    return out


# --------------------------------------------------------------------------
# Settings table  (KEY / TYPE / VALUE)
# --------------------------------------------------------------------------

# LSF:Settings dict key (as lsf_stack.DEFAULT_SETTINGS names it) -> (spec key, TYPE)
_SIMPLE_SETTINGS = {
    "align_transform": ("align.transform", "STR"),
    "align_interp": ("align.interpolation", "STR"),
    "align_minpairs": ("align.min_pairs", "INT"),
    "rejection": ("reject.method", "STR"),
    "weighting": ("weight.scheme", "STR"),
    "norm_method": ("norm.method", "STR"),
    "norm_estimator": ("norm.estimator", "STR"),
    "equalize_rgb": ("norm.equalize_rgb", "BOOL"),
    "skip_normalization": ("norm.skip", "BOOL"),
    "saturation_level": ("norm.saturation_level", "FLT"),
    "output_norm": ("outnorm.requested", "BOOL"),
}
_REJECT_PARAM_KEYS = ("sigma_low", "sigma_high", "max_iters", "gesdt_alpha", "gesdt_max_frac")
_KNOWN_SETTING_KEYS = (
    {v[0] for v in _SIMPLE_SETTINGS.values()}
    | {"sensor.class", "sensor.bayer", "cfa.mode", "grid.oversample", "grid.pixfrac", "grid.kernel",
       "reject.params", "norm.baseline", "outnorm.enabled", "outnorm.min", "outnorm.max",
       "filter", "created_utc", "writer"})


def _fmt_value(typ: str, v) -> str:
    if typ == "BOOL":
        return "T" if v else "F"
    if typ == "INT":
        return str(int(v))
    if typ == "FLT":
        return repr(float(v))
    if typ == "JSON":
        return json.dumps(v, default=str)
    return str(v)


def _parse_value(typ: str, s: str):
    typ = typ.strip().upper()
    s = s if s is not None else ""
    if typ == "BOOL":
        return s.strip().upper() in ("T", "TRUE", "1")
    if typ == "INT":
        return int(float(s))
    if typ == "FLT":
        return float(s)
    if typ == "JSON":
        return json.loads(s) if s.strip() else None
    return s


def _settings_rows(settings: dict, grid: dict, outnorm: dict, baseline, created: str,
                   writer: str, extra_rows: list) -> list:
    """List of (KEY, TYPE, VALUE-string) rows, spec keys first, then unknown/extra ones."""
    rows = {}

    def put(key, typ, val):
        if val is not None:
            rows[key] = (typ, _fmt_value(typ, val))

    sensor_class = grid.get("sensor_class", "mono")
    put("sensor.class", "STR", sensor_class)
    put("sensor.bayer", "STR", grid.get("bayer_pattern") or ("none" if sensor_class == "mono" else "unknown"))
    put("cfa.mode", "STR", _DMOD_FROM_METHOD.get(grid.get("demosaic_method", "mono"), "NONE"))
    put("grid.oversample", "INT", int(grid.get("oversample", 1)))
    put("grid.pixfrac", "FLT", float(grid.get("pixfrac", 1.0)))
    put("grid.kernel", "STR", grid.get("kernel", "none"))
    # Section 6: `filter` is required for mono files (empty string when the frames carried no FILTER).
    put("filter", "STR", grid.get("filter") or ("" if sensor_class == "mono" else None))
    put("reject.params", "JSON", {k: settings.get(k) for k in _REJECT_PARAM_KEYS})
    for k, (key, typ) in _SIMPLE_SETTINGS.items():
        put(key, typ, settings.get(k))
    for k, v in settings.items():             # anything the stacker adds later survives as ext.<name>
        if k not in _SIMPLE_SETTINGS and k not in _REJECT_PARAM_KEYS:
            put("ext." + k, "JSON", v)
    put("outnorm.enabled", "BOOL", bool(outnorm.get("enabled")))
    if outnorm.get("enabled"):
        put("outnorm.min", "FLT", float(outnorm["mini"]))
        put("outnorm.max", "FLT", float(outnorm["maxi"]))
    if baseline:
        put("norm.baseline", "JSON", baseline)
    put("created_utc", "STR", created)
    put("writer", "STR", writer)
    out = [(k, t, v) for k, (t, v) in rows.items()]
    for k, t, v in extra_rows or []:           # C7: unknown rows, unless this write set the same key
        if k not in rows:
            out.append((k, t, v))
    return out


def _props_from_settings(rows: list, phdr, planes: int, shape_hw: tuple) -> tuple:
    """(settings dict, canonical grid dict, outnorm dict, baseline or None, meta extras)."""
    typed, raw = {}, {}
    for key, typ, val in rows:
        raw[key] = (typ, val)
        try:
            typed[key] = _parse_value(typ, val)
        except (ValueError, json.JSONDecodeError):
            typed[key] = val

    settings = {}
    for k, (key, _typ) in _SIMPLE_SETTINGS.items():
        if key in typed:
            settings[k] = typed[key]
    settings.setdefault("saturation_level", None)         # an unset optional is not stored as a row
    for k, v in (typed.get("reject.params") or {}).items():
        settings[k] = v
    for key, v in typed.items():
        if key.startswith("ext."):
            settings[key[4:]] = v

    sens = str(typed.get("sensor.class") or phdr.get("LSFSENS", "mono")).lower()
    bayer = typed.get("sensor.bayer") or (str(phdr.get("LSFBAYR")) if phdr.get("LSFBAYR") else
                                          ("none" if sens == "mono" else "unknown"))
    dmod = str(typed.get("cfa.mode") or phdr.get("LSFDMOD", "NONE")).upper()
    wcs = {k: _py(phdr[k]) for k in ("CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1",
                                      "CD2_2", "CDELT1", "CDELT2", "CROTA2", "CTYPE1", "CTYPE2",
                                      "EQUINOX", "PLTSOLVD") if k in phdr}
    grid = {"width": int(shape_hw[1]), "height": int(shape_hw[0]), "channels": planes, "wcs": wcs,
            "filter": typed.get("filter") or None, "sensor_class": sens, "bayer_pattern": bayer,
            "demosaic_method": _METHOD_FROM_DMOD.get(dmod, "mono"),
            "oversample": typed.get("grid.oversample", 1), "pixfrac": typed.get("grid.pixfrac", 1.0),
            "kernel": typed.get("grid.kernel", "none")}

    # Output normalization lives in the primary header (a merge must not depend on a table).
    norm_on = phdr.get("LSFNORM") is True
    outnorm = {"enabled": False}
    if norm_on and "LSFNMIN" in phdr and "LSFNMAX" in phdr:
        outnorm = {"enabled": True, "mini": float(phdr["LSFNMIN"]), "maxi": float(phdr["LSFNMAX"])}

    extras = [(k, t, v) for k, (t, v) in raw.items() if k not in _KNOWN_SETTING_KEYS
              and not k.startswith("ext.")]
    return settings, grid, outnorm, typed.get("norm.baseline"), extras, typed


# --------------------------------------------------------------------------
# Variable-length string columns (1PA)
#
# astropy hands a PA column back as an array of single characters and drops the spaces,
# which would corrupt file names and verbatim FITS headers. So the heap is parsed here.
# --------------------------------------------------------------------------

def _read_vla_strings(path: str, hdu, name: str) -> list:
    info = hdu.fileinfo()
    nrows, rowlen = int(hdu.header["NAXIS2"]), int(hdu.header["NAXIS1"])
    if nrows == 0:
        return []
    theap = int(hdu.header.get("THEAP", rowlen * nrows))
    q = str(hdu.columns[name].format).upper().split("(")[0].lstrip("0123456789").startswith("Q")
    off = hdu.data.dtype.fields[name][1]
    with open(path, "rb") as f:
        f.seek(info["datLoc"])
        raw = f.read(info["datSpan"])
    wsz = 8 if q else 4
    rows = np.frombuffer(raw, dtype=np.uint8, count=nrows * rowlen).reshape(nrows, rowlen)
    desc = np.ascontiguousarray(rows[:, off:off + 2 * wsz]).view(">i8" if q else ">i4").reshape(nrows, 2)
    return [raw[theap + int(o):theap + int(o) + int(c)].decode("utf-8", "replace") for c, o in desc]


def _vla_column(name: str, values: list):
    arr = np.empty(len(values), dtype=object)
    for i, v in enumerate(values):
        arr[i] = ("" if v is None else str(v)).encode("utf-8", "replace")
    return fits.Column(name=name, format="PA()", array=arr)


def _fixed_str_column(name: str, values: list, min_width: int):
    enc = [("" if v is None else str(v)).encode("ascii", "replace") for v in values]
    w = max([min_width] + [len(e) for e in enc])      # never truncate: widen past the spec minimum
    return fits.Column(name=name, format=f"{w}A", array=np.array(enc, dtype=f"S{w}"))


# --------------------------------------------------------------------------
# Ledger table
# --------------------------------------------------------------------------

_STAT_COLS = (("STAT_MEAN", "mean"), ("STAT_SIG", "sigma"), ("STAT_MAD", "mad"),
              ("STAT_SQBW", "sqrtbwmv"), ("STAT_LOC", "location"), ("STAT_SCL", "scale"))
_KNOWN_LEDGER_COLS = {
    "FRAME_ID", "BATCH", "GEN", "ADDED_UTC", "FILENAME", "DATE_OBS", "EXPTIME", "GAIN", "OFFSET",
    "CCD_TEMP", "FILTER", "ORIGHDR", "REG_TYPE", "REG_XFRM", "ALIGN", "NORM_SC", "NORM_OFF", "NORM_MTH",
    "NORM_EST", "WEIGHT", "FWHM", "WFWHM", "ROUND", "REGQ", "REGBG", "NSTARS", "BKG_MED", "NOISE",
    "NOISE_N"} | {c for c, _ in _STAT_COLS}
_METHOD_W = 12


def _f(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return float("nan")
    return x


def _i(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return INT_NULL
    if not math.isfinite(x) or abs(x) >= 2 ** 31 - 1:
        return INT_NULL
    return int(round(x))


def _plane_list(v, planes: int) -> list:
    if isinstance(v, (list, tuple)):
        out = list(v)
    elif v is None:
        out = []
    else:
        out = [v]
    return (out + [None] * planes)[:planes]


def _float_col(name, vals, fmt="D"):
    dt = ">f8" if fmt == "D" else ">f4"
    return fits.Column(name=name, format=fmt, array=np.array([_f(v) for v in vals], dtype=dt))


def _plane_float_col(name, vals, planes, fmt="D"):
    dt = ">f8" if fmt == "D" else ">f4"
    arr = np.array([[_f(x) for x in _plane_list(v, planes)] for v in vals], dtype=dt)
    return fits.Column(name=name, format=f"{planes}{fmt}", array=arr if planes > 1 else arr[:, 0])


def ledger_to_columns(rows: list, planes: int) -> list:
    """Frame-ledger dicts (the v0.2 field names) -> astropy Columns of LSF_LEDGER."""
    q = lambda r: (r.get("quality") or {})                      # noqa: E731
    nm = lambda r: (r.get("normalization") or {})               # noqa: E731
    al = lambda r: r.get("alignment")                           # noqa: E731
    st = lambda r, k: ((q(r).get("stats") or {}).get(k))        # noqa: E731
    cols = [
        _fixed_str_column("FRAME_ID", [r.get("frame_id") for r in rows], 32),
        fits.Column(name="BATCH", format="J", null=INT_NULL,
                    array=np.array([_i(r.get("batch")) for r in rows], dtype=">i4")),
        fits.Column(name="GEN", format="J", null=INT_NULL,
                    array=np.array([_i(r.get("gen")) for r in rows], dtype=">i4")),
        _fixed_str_column("ADDED_UTC", [r.get("added_utc") for r in rows], 24),
        _vla_column("FILENAME", [r.get("original_filename") for r in rows]),
        _fixed_str_column("DATE_OBS", [r.get("capture_timestamp") for r in rows], 24),
        _float_col("EXPTIME", [r.get("exposure_seconds") for r in rows], "D"),
        fits.Column(name="GAIN", format="J", null=INT_NULL,
                    array=np.array([_i(r.get("gain")) for r in rows], dtype=">i4")),
        fits.Column(name="OFFSET", format="J", null=INT_NULL,
                    array=np.array([_i(r.get("offset")) for r in rows], dtype=">i4")),
        _float_col("CCD_TEMP", [r.get("sensor_temp_c") for r in rows], "E"),
        _fixed_str_column("FILTER", [r.get("filter") for r in rows], 16),
        _vla_column("ORIGHDR", [r.get("original_header") for r in rows]),
        _fixed_str_column("REG_TYPE", [(al(r) or {}).get("transform") if isinstance(al(r), dict) else None
                                        for r in rows], 8),
    ]
    xf = np.full((len(rows), 3, 3), np.nan, dtype=">f8")
    for n, r in enumerate(rows):
        h = r.get("registration_transform")
        if isinstance(h, (list, tuple)) and len(h) == 9:
            xf[n] = np.array([_f(x) for x in h], dtype=float).reshape(3, 3)
    cols.append(fits.Column(name="REG_XFRM", format="9D", dim="(3,3)", array=xf))
    cols.append(_vla_column("ALIGN", [json.dumps(al(r), default=str) if al(r) is not None else "" for r in rows]))
    cols.append(_plane_float_col("NORM_SC", [nm(r).get("scale") for r in rows], planes, "D"))
    cols.append(_plane_float_col("NORM_OFF", [nm(r).get("offset") for r in rows], planes, "D"))
    meth = np.array([[("" if m is None else str(m)).encode("ascii", "replace")[:_METHOD_W]
                      for m in _plane_list(nm(r).get("method"), planes)] for r in rows], dtype=f"S{_METHOD_W}")
    if planes > 1:
        cols.append(fits.Column(name="NORM_MTH", format=f"{planes * _METHOD_W}A",
                                dim=f"({_METHOD_W},{planes})", array=meth))
    else:
        cols.append(fits.Column(name="NORM_MTH", format=f"{_METHOD_W}A", array=meth[:, 0]))
    # Not in the spec's column table: the addscale estimator (kind, loc, spread) of v0.2 Section 8.
    cols.append(_vla_column("NORM_EST", [json.dumps(nm(r)["estimator"], default=str)
                                         if nm(r).get("estimator") is not None else "" for r in rows]))
    cols.append(_plane_float_col("WEIGHT", [r.get("weight_applied") for r in rows], planes, "D"))
    cols += [_float_col("FWHM", [q(r).get("fwhm") for r in rows], "E"),
             _float_col("WFWHM", [q(r).get("weighted_fwhm") for r in rows], "E"),
             _float_col("ROUND", [q(r).get("roundness") for r in rows], "E"),
             _float_col("REGQ", [q(r).get("registration_quality") for r in rows], "E"),
             _float_col("REGBG", [q(r).get("registration_background") for r in rows], "E"),
             fits.Column(name="NSTARS", format="J", null=INT_NULL,
                         array=np.array([_i(q(r).get("star_count")) for r in rows], dtype=">i4")),
             _plane_float_col("BKG_MED", [q(r).get("background_median") for r in rows], planes, "E"),
             _plane_float_col("NOISE", [q(r).get("noise") for r in rows], planes, "E"),
             _plane_float_col("NOISE_N", [q(r).get("noise_normalized") for r in rows], planes, "E")]
    for col, key in _STAT_COLS:
        cols.append(_plane_float_col(col, [st(r, key) for r in rows], planes, "E"))
    return cols


def _nan_none(x):
    x = float(x)
    return None if math.isnan(x) else x


def ledger_from_hdu(path: str, hdu, planes: int):
    """LSF_LEDGER -> (list of ledger dicts, list of unknown-column records)."""
    data = hdu.data
    names = {c.name.upper(): c.name for c in hdu.columns}
    n = len(data)

    def col(name):
        return data[names[name]] if name in names else None

    def vla(name):
        return _read_vla_strings(path, hdu, names[name]) if name in names else [None] * n

    def strs(name):
        c = col(name)
        return [None] * n if c is None else [(str(x).strip() or None) for x in c]

    def ints(name):
        c = col(name)
        if c is None:
            return [None] * n
        null = hdu.columns[names[name]].null
        null = INT_NULL if null is None else null
        return [None if int(x) == null else int(x) for x in c]

    def flts(name):
        c = col(name)
        return [None] * n if c is None else [_nan_none(x) for x in c]

    def plane_lists(name):
        c = col(name)
        if c is None:
            return [None] * n
        return [[_nan_none(x) for x in np.atleast_1d(row)] for row in c]

    frame_id, added, date_obs, filt = strs("FRAME_ID"), strs("ADDED_UTC"), strs("DATE_OBS"), strs("FILTER")
    batch, gen, gain, offset, nstars = ints("BATCH"), ints("GEN"), ints("GAIN"), ints("OFFSET"), ints("NSTARS")
    fname, orighdr, align_s, est_s = vla("FILENAME"), vla("ORIGHDR"), vla("ALIGN"), vla("NORM_EST")
    exptime, ccd = flts("EXPTIME"), flts("CCD_TEMP")
    fwhm, wfwhm, rnd, regq, regbg = flts("FWHM"), flts("WFWHM"), flts("ROUND"), flts("REGQ"), flts("REGBG")
    reg_type = strs("REG_TYPE")
    xfrm = col("REG_XFRM")
    sc, off, wt = plane_lists("NORM_SC"), plane_lists("NORM_OFF"), plane_lists("WEIGHT")
    bkg, noise, noise_n = plane_lists("BKG_MED"), plane_lists("NOISE"), plane_lists("NOISE_N")
    stats = {key: plane_lists(c) for c, key in _STAT_COLS}
    mc = col("NORM_MTH")
    methods = [None] * n if mc is None else [[(str(x).strip() or None) for x in np.atleast_1d(r)] for r in mc]

    rows = []
    for k in range(n):
        h = None
        if xfrm is not None:
            flat = np.asarray(xfrm[k], dtype=float).ravel()
            h = None if np.isnan(flat).any() else [float(x) for x in flat]
        alignment = None
        if align_s[k]:
            try:
                alignment = json.loads(align_s[k])
            except json.JSONDecodeError:
                alignment = None
        if alignment is None and reg_type[k]:
            alignment = {"transform": reg_type[k]}
        norm = {"scale": sc[k], "offset": off[k], "method": methods[k]}
        if est_s[k]:
            try:
                norm["estimator"] = json.loads(est_s[k])
            except json.JSONDecodeError:
                pass
        rows.append({
            "frame_id": frame_id[k], "batch": batch[k], "gen": gen[k], "added_utc": added[k],
            "original_filename": fname[k] or None, "capture_timestamp": date_obs[k],
            "exposure_seconds": exptime[k], "gain": gain[k], "offset": offset[k],
            "sensor_temp_c": ccd[k], "filter": filt[k], "original_header": orighdr[k] or None,
            "registration_transform": h, "alignment": alignment, "normalization": norm,
            "weight_applied": wt[k],
            "quality": {"fwhm": fwhm[k], "weighted_fwhm": wfwhm[k], "roundness": rnd[k],
                        "star_count": nstars[k], "registration_quality": regq[k],
                        "registration_background": regbg[k], "background_median": bkg[k],
                        "noise": noise[k], "noise_normalized": noise_n[k],
                        "stats": {key: stats[key][k] for key in stats}},
        })

    unknown = []                         # C7: columns this reader does not know
    for c in hdu.columns:
        if c.name.upper() in _KNOWN_LEDGER_COLS:
            continue
        fmt = str(c.format)
        base = fmt.upper().split("(")[0].lstrip("0123456789")
        if base[:1] in ("P", "Q"):
            if base[1:2] == "A":
                arr = np.array([t.encode("utf-8") for t in _read_vla_strings(path, hdu, c.name)], dtype=object)
                vla_kind = "str"
            else:
                arr = [np.array(x) for x in data[c.name]]
                vla_kind = True
        else:
            arr, vla_kind = np.array(data[c.name]), False
        unknown.append({"name": c.name, "format": fmt, "dim": c.dim, "null": c.null, "unit": c.unit,
                        "array": arr, "vla": vla_kind})
    return rows, unknown


def _pad_unknown(rec: dict, n_old: int, n_total: int):
    """Rebuild a preserved unknown column for a ledger that gained rows: new rows get the
    spec's 'unknown' value (NaN for floats, TNULL for integers, empty for strings/arrays)."""
    add = n_total - n_old
    arr, fmt = rec["array"], rec["format"]
    if rec["vla"]:
        is_str = rec["vla"] == "str"
        base = None if is_str else (arr[0].dtype if len(arr) else np.dtype("f8"))
        out = np.empty(n_total, dtype=object)
        for i in range(n_old):
            out[i] = arr[i]
        for i in range(n_old, n_total):
            out[i] = b"" if is_str else np.array([], dtype=base)
        return fits.Column(name=rec["name"], format=re.sub(r"\(\d*\)$", "()", fmt), array=out,
                           unit=rec["unit"])
    if arr.dtype.kind == "f":
        fill = np.full((add,) + arr.shape[1:], np.nan, dtype=arr.dtype)
    elif arr.dtype.kind in "iu":
        nv = rec["null"] if rec["null"] is not None else 0
        fill = np.full((add,) + arr.shape[1:], nv, dtype=arr.dtype)
    elif arr.dtype.kind in "SU":
        fill = np.full((add,) + arr.shape[1:], "", dtype=arr.dtype)
    else:
        fill = np.zeros((add,) + arr.shape[1:], dtype=arr.dtype)
    return fits.Column(name=rec["name"], format=fmt, dim=rec["dim"], unit=rec["unit"],
                       null=rec["null"], array=np.concatenate([arr, fill]) if add else arr)


# --------------------------------------------------------------------------
# Reader
# --------------------------------------------------------------------------

def _hdu_data(path: str, hdu):
    """hdu.data, memory-mapped when the HDU is unscaled. An HDU with BZERO/BSCALE (the uint16
    count plane, or a primary written by a tool that adds BSCALE = 1 / BZERO = 0) cannot be
    mapped by astropy; it is read into RAM instead (the count plane is 2 bytes per pixel)."""
    try:
        return hdu.data
    except ValueError:
        idx = None
        with fits.open(path, memmap=False) as h2:
            for i, h in enumerate(h2):
                if i and _hdu_name(h) == _hdu_name(hdu):
                    idx = i
                    break
            if idx is None:
                return np.array(fits.open(path, memmap=False)[0].data)
            return np.array(h2[idx].data)


def read_lsf_fits(path: str) -> dict:
    """Open an LSF FITS file and return:

        {"images": {"Stack"|"WeightSum"|"SumSq"|"SumW2"|"ContributionCount":
                        {"data": channel-last array, "fits_keywords": [(key, value, comment)]}},
         "properties": {"LSF:FrameLedger", "LSF:CanonicalGrid", "LSF:Settings", "LSF:OutputNorm",
                        ["LSF:NormBaseline"]}   # JSON strings
         "meta": {...}, "format": "fits"}

    Uncompressed float planes are read-only memory maps (big-endian; use native() before
    handing a slice to torch or any code that needs native byte order). Compressed planes
    are decompressed into RAM. No validation is done here (a foreign file must stay
    readable, C1): call validate_lsf_fits() before merging.
    """
    _require_astropy()
    hdul = fits.open(path, memmap=True)
    try:
        prim = hdul[0]
        phdr = prim.header
        pdata = _hdu_data(path, prim) if not (prim.header.get("BSCALE") or prim.header.get("BZERO")) \
            else np.array(fits.getdata(path, 0, memmap=False))
        if pdata is None or pdata.ndim not in (2, 3):
            raise LsfFitsError(f"'{path}' has no 2-D / 3-D image in HDU 0.")
        exts = _find_hdus(hdul)
        # LSF_W2 is the one plane that may be absent: a file written before it existed is still a
        # readable image (C1) and its other planes can still be exported. It cannot be merged into
        # (validate_hdul check 3), and lsf_common.check_lsf_file() flags the missing 'SumW2' image.
        missing = [n for n in KNOWN_EXTS if n not in exts and n != EXT_W2]
        if missing:
            raise LsfFitsError(f"'{path}' is not a valid LSF file (missing HDU: {', '.join(missing)}).")
        missing_hdus = [n for n in KNOWN_EXTS if n not in exts]

        stack = _from_fits_planes(pdata)
        planes = stack.shape[2] if stack.ndim == 3 else 1
        images = {"Stack": {"data": stack, "fits_keywords": [
            (c.keyword, _py(c.value), c.comment) for c in phdr.cards
            if not _STRUCTURAL_RE.match(c.keyword.upper())]}}
        for mem, name in IMAGE_EXTS.items():
            if name not in exts:                       # only LSF_W2 can get here
                continue
            arr = _from_fits_planes(_hdu_data(path, exts[name]))
            if mem == "ContributionCount" and stack.ndim == 3 and arr.ndim == 2:
                arr = np.broadcast_to(arr[..., None], stack.shape)       # shared count (NAXIS3 = 1)
            elif mem == "ContributionCount" and arr.ndim == 3 and arr.shape[2] == 1 and stack.ndim == 3:
                arr = np.broadcast_to(arr, stack.shape)
            images[mem] = {"data": arr, "fits_keywords": []}

        # settings
        sh = exts[EXT_SETTINGS]
        keys = [str(x).strip() for x in sh.data["KEY"]]
        types = [str(x).strip() for x in sh.data["TYPE"]]
        vals = _read_vla_strings(path, sh, "VALUE")
        settings, grid, outnorm, baseline, s_extra, typed = _props_from_settings(
            list(zip(keys, types, vals)), phdr, planes, stack.shape[:2])

        ledger, unknown_cols = ledger_from_hdu(path, exts[EXT_LEDGER], planes)

        props = {"LSF:FrameLedger": json.dumps(ledger, default=str, indent=2),
                 "LSF:CanonicalGrid": json.dumps(grid, default=str, indent=2),
                 "LSF:Settings": json.dumps(settings, default=str, indent=2),
                 "LSF:OutputNorm": json.dumps(outnorm, default=str)}
        if baseline:
            props["LSF:NormBaseline"] = json.dumps(baseline, default=str, indent=2)

        # C7: everything a rewrite must carry forward. Unknown HDUs are kept as raw bytes
        # (header text + padded data), so even compressed or heap-carrying foreign HDUs survive exactly.
        extra_hdus = []
        for i, h in enumerate(hdul):
            if i == 0 or _hdu_name(h) in KNOWN_EXTS:
                continue
            info = h.fileinfo()
            with open(path, "rb") as fh:
                fh.seek(info["hdrLoc"])
                head = fh.read(info["datLoc"] - info["hdrLoc"])
                fh.seek(info["datLoc"])
                body = fh.read(info["datSpan"])
            extra_hdus.append({"name": _hdu_name(h), "header": head, "data": body})
        meta = {
            "uid": phdr.get("LSFUID"), "gen": int(phdr.get("LSFGEN", 0) or 0),
            "version": phdr.get("LSFVER"), "row_order": phdr.get("ROWORDER"),
            "created_utc": typed.get("created_utc"), "writer": typed.get("writer"),
            "primary_cards": [c.image for c in phdr.cards if not _STRUCTURAL_RE.match(c.keyword.upper())],
            "settings_extra": s_extra, "extra_hdus": extra_hdus,
            "ledger_unknown_columns": unknown_cols, "n_ledger_rows": len(ledger),
            "missing_hdus": missing_hdus,
        }
        return {"images": images, "properties": props, "meta": meta, "format": "fits"}
    finally:
        hdul.close()          # the memmaps stay valid while referenced; drop them before overwriting


# --------------------------------------------------------------------------
# Validation (Section 8)
# --------------------------------------------------------------------------

def _close(a, b, rel=1e-12) -> bool:
    return math.isclose(float(a), float(b), rel_tol=rel, abs_tol=0.0)


def validate_hdul(hdul) -> None:
    """Run the six Section 8 pre-merge checks, in order. Raises LsfValidationError naming
    the first one that fails. Reads the primary data once (check 5)."""
    phdr = hdul[0].header

    # 1. identity and version
    if phdr.get("LSF") is not True:
        raise LsfValidationError(1, "the primary header has no 'LSF = T' (not an LSF file, or the "
                                    "LSF keywords were stripped)")
    ver = str(phdr.get("LSFVER", "")).strip()
    try:
        major = int(ver.split(".")[0])
    except ValueError:
        raise LsfValidationError(1, f"LSFVER '{ver}' is not readable") from None
    if major != SUPPORTED_MAJOR:
        raise LsfValidationError(1, f"LSFVER {ver}: major version {major} is not supported "
                                    f"(this writer handles {SUPPORTED_MAJOR}.x)")

    # 2. interrupted update
    if phdr.get("LSFBUSY") is True:
        raise LsfValidationError(2, "LSFBUSY = T: an in-place update was interrupted; restore a backup "
                                    "(the Stack is still viewable)")

    # 3. required HDUs, found by name, with consistent dimensions
    exts = _find_hdus(hdul)
    missing = [n for n in KNOWN_EXTS if n not in exts]
    if missing:
        hint = ("a tool that saved only the primary, or dropped the extensions"
                + ("; LSF_W2 is the sum-of-squared-weights plane added to the profile in v0.3, so a file "
                   "written before it existed cannot be merged into and the stack must be rebuilt from "
                   "the original frames" if EXT_W2 in missing else ""))
        raise LsfValidationError(3, f"missing HDU(s) {', '.join(missing)} ({hint})")
    pshape = tuple(hdul[0].shape)
    for n in (EXT_WEIGHT, EXT_M2, EXT_W2):           # all per channel, same grid as the Stack
        if tuple(exts[n].shape) != pshape:
            raise LsfValidationError(3, f"{n} has shape {tuple(exts[n].shape)}, the Stack has {pshape}")
    cshape = tuple(exts[EXT_COUNT].shape)
    if cshape != pshape and not (len(pshape) == 3 and cshape == (1,) + pshape[1:]):
        raise LsfValidationError(3, f"LSF_COUNT has shape {cshape}, the Stack has {pshape}")
    ro = phdr.get("ROWORDER")
    for n in (EXT_WEIGHT, EXT_M2, EXT_W2, EXT_COUNT):
        r = exts[n].header.get("ROWORDER")
        if r is not None and ro is not None and str(r) != str(ro):
            raise LsfValidationError(3, f"{n} ROWORDER '{r}' differs from the primary's '{ro}'")

    # 4. binding: LSFUID / LSFGEN identical in every LSF HDU
    uid, gen = phdr.get("LSFUID"), phdr.get("LSFGEN")
    if uid is None or gen is None:
        raise LsfValidationError(4, "the primary header lacks LSFUID or LSFGEN")
    for name, h in exts.items():
        if not name.startswith("LSF_"):
            continue
        if h.header.get("LSFUID") != uid or h.header.get("LSFGEN") != gen:
            raise LsfValidationError(
                4, f"{name} has LSFUID/LSFGEN {h.header.get('LSFUID')}/{h.header.get('LSFGEN')}, the primary "
                   f"has {uid}/{gen} (the file was rewritten by a tool that kept only part of it)")

    # 5. the primary's data is the data that was committed
    want = phdr.get("LSFDSUM")
    if want is None:
        raise LsfValidationError(5, "the primary header lacks LSFDSUM")
    data = hdul[0].data
    have = datasum_of_planes(data) if data.dtype.itemsize == 4 else None
    if have is None or str(have) != str(want).strip():
        raise LsfValidationError(5, f"DATASUM of the primary data is {have}, LSFDSUM says {want} "
                                    f"(the Stack was changed after the last LSF commit)")

    # 6. settings mirrored in the primary header agree with the table
    sh = exts[EXT_SETTINGS]
    typed = {}
    keys = [str(x).strip() for x in sh.data["KEY"]]
    types = [str(x).strip() for x in sh.data["TYPE"]]
    vals = _read_vla_strings(hdul.filename(), sh, "VALUE") if hdul.filename() else [""] * len(keys)
    for k, t, v in zip(keys, types, vals):
        typed[k] = _parse_value(t, v)

    def bad(what, hv, tv):
        raise LsfValidationError(6, f"{what}: primary header says {hv!r}, LSF_SETTINGS says {tv!r}")

    if "sensor.class" in typed and str(typed["sensor.class"]).upper() != str(phdr.get("LSFSENS", "")).upper():
        bad("sensor class (LSFSENS / sensor.class)", phdr.get("LSFSENS"), typed["sensor.class"])
    if "LSFBAYR" in phdr and "sensor.bayer" in typed \
            and str(typed["sensor.bayer"]).upper() != str(phdr["LSFBAYR"]).upper():
        bad("Bayer pattern (LSFBAYR / sensor.bayer)", phdr["LSFBAYR"], typed["sensor.bayer"])
    if "cfa.mode" in typed and str(typed["cfa.mode"]).upper() != str(phdr.get("LSFDMOD", "")).upper():
        bad("demosaic mode (LSFDMOD / cfa.mode)", phdr.get("LSFDMOD"), typed["cfa.mode"])
    if "outnorm.enabled" in typed and bool(typed["outnorm.enabled"]) != (phdr.get("LSFNORM") is True):
        bad("output normalization (LSFNORM / outnorm.enabled)", phdr.get("LSFNORM"), typed["outnorm.enabled"])
    if phdr.get("LSFNORM") is True:
        for hk, sk in (("LSFNMIN", "outnorm.min"), ("LSFNMAX", "outnorm.max")):
            if hk not in phdr:
                raise LsfValidationError(6, f"{hk} is missing from a Stack flagged LSFNORM = T; the "
                                            f"normalization cannot be inverted")
            if sk not in typed or not _close(phdr[hk], typed[sk]):
                bad(f"{hk} / {sk}", phdr[hk], typed.get(sk))


def validate_lsf_fits(path: str) -> None:
    """Section 8 pre-merge validation of the file at `path`; raises LsfValidationError."""
    _require_astropy()
    with fits.open(path, memmap=True) as hdul:
        validate_hdul(hdul)


# --------------------------------------------------------------------------
# Writer
# --------------------------------------------------------------------------

def _wcs_cards(phdr) -> list:
    return [c.image for c in phdr.cards if is_wcs_key(c.keyword)]


def _image_hdu(name, arr, uid, gen, row_order, wcs_images, compress, extra=None):
    data = _to_fits_planes(arr)
    hdr = fits.Header()
    hdr["EXTVER"] = (1, "LSF extension version")
    hdr["ROWORDER"] = (row_order, "same row order as the primary")
    hdr["LSFUID"] = (uid, "LSF file id")
    hdr["LSFGEN"] = (gen, "LSF generation counter")
    for k, (v, c) in (extra or {}).items():
        hdr[k] = (v, c)
    for img in wcs_images:                    # each image carries its own WCS (INHERIT is unevenly supported)
        hdr.append(fits.Card.fromstring(img), end=True)
    if compress:
        ctype = "RICE_1" if data.dtype.kind == "u" else "GZIP_2"
        kw = {} if data.dtype.kind == "u" else {"quantize_level": 0}      # lossless: no float quantization
        hdu = fits.CompImageHDU(data=data, header=hdr, name=name, compression_type=ctype, **kw)
    else:
        hdu = fits.ImageHDU(data=data, header=hdr, name=name, ver=1)
    return hdu


def _append_hdu(path: str, hdu, post_checksum: bool) -> None:
    """Append one HDU. Plain images can be checksummed before they are written. Tables with a
    heap (the variable-length columns) and compressed images only get a valid CHECKSUM/DATASUM
    once astropy has serialized them, so those are checksummed in a second pass over the HDU
    that was just appended (always the last one)."""
    if not post_checksum:
        hdu.add_checksum()
    with fits.open(path, mode="append") as hl:
        hl.append(hdu)
    if post_checksum:
        with fits.open(path, mode="update") as hl:
            hl[-1].add_checksum()


def write_lsf_fits(path: str, images: list, properties: dict, compress=None, meta: dict | None = None,
                   new_rows: int | None = None, row_order: str | None = None,
                   writer: str | None = None) -> dict:
    """Write a complete LSF FITS file at `path` (the caller passes a temporary name and renames
    it over the original afterwards: Section 10's commit protocol).

    images       objects with .id ("Stack", "WeightSum", "SumSq", "SumW2", "ContributionCount"), .data
                 (channel-last (H, W[, C]) arrays) and .fits_keywords [(key, value, comment)] on Stack
    properties   the JSON blocks (LSF:FrameLedger / CanonicalGrid / Settings /
                 OutputNorm / [NormBaseline])
    compress     falsy | True | "gzip2": extension images tile-compressed, losslessly
    meta         the "meta" of read_lsf_fits() for the file being appended to (None on initialize).
                 Carries LSFUID, LSFGEN (+1), created_utc and every unknown thing to preserve.
    new_rows     how many trailing ledger rows are new in this commit (stamped with this commit's
                 generation / batch / time). Rows before them that carry no 'gen' are legacy (0).
    """
    _require_astropy()
    meta = meta or {}
    imgs = {im.id: im for im in images}
    for need in ("Stack",) + tuple(IMAGE_EXTS):
        if need not in imgs:
            raise LsfFitsError(f"write_lsf_fits: image '{need}' is missing")
    stack = imgs["Stack"].data
    if not (stack.dtype.kind == "f" and stack.dtype.itemsize == 4):      # big-endian memmaps are fine
        raise LsfFitsError("The Stack must be float32 (C1)")
    if stack.ndim == 3 and stack.shape[2] != 3:
        raise LsfFitsError(f"Stack has {stack.shape[2]} channels; the profile allows mono or 3-plane RGB")
    planes = stack.shape[2] if stack.ndim == 3 else 1

    ledger = json.loads(properties["LSF:FrameLedger"])
    grid = json.loads(properties["LSF:CanonicalGrid"])
    settings = json.loads(properties["LSF:Settings"])
    outnorm = json.loads(properties["LSF:OutputNorm"])
    baseline = json.loads(properties["LSF:NormBaseline"]) if properties.get("LSF:NormBaseline") else None

    uid = meta.get("uid") or uuid.uuid4().hex
    gen = int(meta.get("gen") or 0) + 1
    created = meta.get("created_utc") or _now_utc()
    now = _now_utc()
    row_order = row_order or meta.get("row_order") or DEFAULT_ROW_ORDER
    writer = writer or f"lsf_fits.py {PROFILE_VERSION}"
    compress = bool(compress)

    # ---- ledger bookkeeping: gen / batch / added_utc ----
    n_new = len(ledger) if new_rows is None else int(new_rows)
    n_old = len(ledger) - n_new
    batch_no = max([int(r["batch"]) for r in ledger[:n_old] if r.get("batch") is not None] or [0]) + 1
    for i, r in enumerate(ledger):
        if i >= n_old:
            r["gen"], r["batch"], r["added_utc"] = gen, batch_no, now
        else:
            r.setdefault("gen", 0)
            r.setdefault("batch", 0)
            r.setdefault("added_utc", None)

    # ---- primary header ----
    sens = "CFA" if grid.get("sensor_class") == "cfa" else "MONO"
    hdr = fits.Header()
    hdr["ROWORDER"] = (row_order, "row order of the pixel array")
    provided = {}
    commentary_new = []
    for key, value, comment in imgs["Stack"].fits_keywords:
        k = str(key).upper()
        if _STRUCTURAL_RE.match(k) or k in LSF_ID_KEYS or k in BAYER_KEYS or value is None:
            continue
        if k in _COMMENTARY:
            commentary_new.append((k, str(value)))
        else:
            provided[k] = (_py(value), comment or "")
    old_cards = [fits.Card.fromstring(img) for img in meta.get("primary_cards", [])]
    for c in old_cards:                                       # C7: unknown keywords survive
        k = c.keyword.upper()
        if k in _COMMENTARY or k in provided or k in LSF_ID_KEYS or k in BAYER_KEYS or k in ("LSFCSAT",):
            continue
        hdr.append(c, end=True)
    for k, (v, c) in provided.items():
        try:
            hdr[k] = (v, c)
        except (ValueError, TypeError):
            pass
    hdr["LSF"] = (True, "Living Stack Format file")
    hdr["LSFVER"] = (PROFILE_VERSION, "LSF FITS profile version")
    hdr["LSFUID"] = (uid, "unique id, never changes")
    hdr["LSFGEN"] = (gen, "generation counter, +1 per commit")
    hdr["LSFSENS"] = (sens, "sensor class (permanent)")
    bayer = str(grid.get("bayer_pattern") or "").upper()
    if sens == "CFA" and bayer in _BAYER_PATTERNS:
        hdr["LSFBAYR"] = (bayer, "Bayer pattern of the source frames")
    hdr["LSFDMOD"] = (_DMOD_FROM_METHOD.get(grid.get("demosaic_method", "mono"), "NONE"),
                      "demosaic mode (permanent)")
    on = bool(outnorm.get("enabled"))
    hdr["LSFNORM"] = (on, "Stack is output-normalized")
    if on:
        _set_card(hdr, "LSFNMIN", float(outnorm["mini"]), "raw value mapped to 0")
        _set_card(hdr, "LSFNMAX", float(outnorm["maxi"]), "raw value mapped to 1")
    dsum = datasum_of_planes(stack)
    hdr["LSFDSUM"] = (str(dsum), "DATASUM of the primary data at commit")
    for c in old_cards:
        if c.keyword.upper() in _COMMENTARY:
            hdr.append(c, end=True)
    for k, text in commentary_new:
        (hdr.add_history if k == "HISTORY" else hdr.add_comment)(text)
    # the WCS goes into every extension image too
    wcs_images = _wcs_cards(hdr)

    pdata = _to_fits_planes(stack)
    phdu = fits.PrimaryHDU(data=pdata, header=hdr)
    fits.HDUList([phdu]).writeto(path, overwrite=True, checksum=True)
    del phdu, pdata

    # ---- extension images ----
    count = imgs["ContributionCount"].data
    csat = bool(count.size and int(count.max()) >= 65535)
    for mem, name in IMAGE_EXTS.items():
        arr = imgs[mem].data
        want = np.uint16 if mem == "ContributionCount" else np.float32
        if arr.dtype.kind != np.dtype(want).kind or arr.dtype.itemsize != np.dtype(want).itemsize:
            arr = np.asarray(arr, dtype=want)
        extra = {"LSFCSAT": (True, "count saturated at 65535 somewhere")} if (mem == "ContributionCount" and csat) else None
        hdu = _image_hdu(name, arr, uid, gen, row_order, wcs_images, compress, extra)
        _append_hdu(path, hdu, compress)
        del hdu

    # ---- settings table ----
    rows = _settings_rows(settings, grid, outnorm, baseline, created, writer, meta.get("settings_extra"))
    cols = [_fixed_str_column("KEY", [r[0] for r in rows], 32),
            _fixed_str_column("TYPE", [r[1] for r in rows], 8),
            _vla_column("VALUE", [r[2] for r in rows])]
    sh = fits.BinTableHDU.from_columns(cols, name=EXT_SETTINGS)
    sh.header["EXTVER"] = 1
    sh.header["LSFUID"] = (uid, "LSF file id")
    sh.header["LSFGEN"] = (gen, "LSF generation counter")
    _append_hdu(path, sh, True)

    # ---- unknown HDUs from the file being rewritten (C7), appended verbatim.
    # An unknown LSF_ HDU gets the new binding; its CHECKSUM is dropped (its header changed).
    for rec in meta.get("extra_hdus", []):
        head = rec["header"]
        if rec["name"].startswith("LSF_"):
            h = fits.Header.fromstring(head.decode("ascii", "replace"))
            h["LSFUID"], h["LSFGEN"] = uid, gen
            for k in ("CHECKSUM", "DATASUM"):
                if k in h:
                    del h[k]
            head = h.tostring(sep="", endcard=True, padding=True).encode("ascii")
        with open(path, "ab") as f:
            f.write(head)
            f.write(rec["data"])

    # ---- ledger (always last) ----
    lcols = ledger_to_columns(ledger, planes)
    for rec in meta.get("ledger_unknown_columns", []):
        if len(rec["array"]) == len(ledger):
            lcols.append(fits.Column(name=rec["name"], format=rec["format"], dim=rec["dim"],
                                     unit=rec["unit"], null=rec["null"], array=rec["array"]))
        else:
            lcols.append(_pad_unknown(rec, len(rec["array"]), len(ledger)))
    lh = fits.BinTableHDU.from_columns(lcols, name=EXT_LEDGER)
    lh.header["EXTVER"] = 1
    lh.header["LSFUID"] = (uid, "LSF file id")
    lh.header["LSFGEN"] = (gen, "LSF generation counter")
    _append_hdu(path, lh, True)

    # Section 10 step 4: flush to disk before the caller's atomic rename
    with open(path, "rb+") as f:
        f.flush()
        os.fsync(f.fileno())
    return {"bytes": os.path.getsize(path), "mode": "gzip2" if compress else None,
            "uid": uid, "gen": gen}