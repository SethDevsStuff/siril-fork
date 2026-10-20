"""
lsf_common.py -- shared helpers for the Living Stack Format (LSF) proof of concept.

LSF files are FITS files (live_stack_fits.md); the container itself (reader, writer, Section 8
validation) lives in lsf_fits.py. This module holds the container-independent pieces: the
in-memory image type, JSON helpers, the processing/statistics code the stacker and exporter
share, and thin wrappers (read_lsf / write_lsf / validate_lsf_file) over lsf_fits.

This module has no dependency on sirilpy so it can be unit tested on its
own (see test_lsf_fits.py).
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from typing import Optional
from concurrent.futures import ThreadPoolExecutor
import os

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

try:
    from scipy.stats import t as t_dist
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

@dataclass
class LsfImage:
    """One image plane on its way into an LSF file."""
    id: str
    data: np.ndarray  # shape (H, W) or (H, W, C), channel-last in memory
    fits_keywords: list = field(default_factory=list)  # list of (key, value, comment)


def write_fits(path: str, data: np.ndarray, header_pairs: list) -> None:
    """Write a FITS file directly using astropy, bypassing sirilpy.
    
    Args:
        path: Output FITS file path
        data: Image data array (H, W) or (H, W, C)
        header_pairs: List of (key, value, comment) tuples for FITS header
    """
    try:
        from astropy.io import fits
    except ImportError:
        raise RuntimeError("astropy is required for direct FITS writing. "
                          "Install with: pip install astropy")
    
    # Ensure contiguous array for FITS writing
    data = np.ascontiguousarray(data, dtype=np.float32)
    
    # Handle multi-channel data: FITS expects (C, H, W) for multi-plane images
    if data.ndim == 3:
        data = np.transpose(data, (2, 0, 1))
    
    hdu = fits.PrimaryHDU(data)
    
    # Add FITS keywords from header_pairs
    for key, value, comment in header_pairs:
        if key == "HISTORY":
            hdu.header.add_history(str(value))
        elif key == "COMMENT":
            hdu.header.add_comment(str(value))
        else:
            try:
                hdu.header[key] = (value, comment if comment else "")
            except Exception:
                pass
    
    hdu.writeto(path, overwrite=True)


def fits_header_string(cards) -> str:
    """Well-formed FITS header string (80-column cards, END) for
    SirilInterface.save_image_file(). `cards`: (key, value) or (key, value, comment) tuples;
    HISTORY / COMMENT keys become text cards (split into 72-character pieces)."""
    out = []
    for card in cards:
        key, val = card[0], card[1]
        comment = card[2] if len(card) > 2 else ""
        if key in ("HISTORY", "COMMENT"):
            text = str(val)
            for i in range(0, max(len(text), 1), 72):
                out.append(f"{key:<8}{text[i:i + 72]}".ljust(80))
            continue
        if isinstance(val, bool):
            lit = ("T" if val else "F").rjust(20)
        elif isinstance(val, str):
            lit = ("'" + val.replace("'", "''") + "'").ljust(20)
        else:
            lit = str(val).rjust(20)
        text = f"{key:<8}= {lit}"
        if comment:
            text += f" / {comment}"
        out.append(text.ljust(80)[:80])
    out.append("END".ljust(80))
    return "\n".join(out)


# --------------------------------------------------------------------------
# Ingestion-time normalization (live_stack.md v0.2, Section 8.3).
#
# This is a *pixel-pairwise* robust linear regression of an incoming
# frame's values against a reference image's values at the same spatial
# positions ("reference ~= source*scale + offset"), analogous to Siril's
# -norm=addscale or PixInsight's LinearFit run between two *images* --
# fundamentally different from matching each image's own marginal
# location/scale statistics (what this file's earlier robust_location_scale
# + a per-image normalize_frame() used to do for LSF, before v0.2).
#
# That marginal-statistics approach breaks down exactly where LSF needs
# normalization most: once the reference is a many-frame accumulator mean,
# its pixel-to-pixel *spread* (MAD) shrinks roughly like 1/sqrt(N) relative
# to a single incoming frame's noise, even though the underlying flux scale
# hasn't changed at all. Matching marginal spreads in that situation scales
# a real single frame's noise down by that same factor -- silently
# corrupting it. A pixel-pairwise fit has no such failure mode: it compares
# the two images point-by-point, not their independent marginal statistics.
# --------------------------------------------------------------------------

def valid_pixel_mask(frame: np.ndarray, saturation_level: "float | None" = None) -> np.ndarray:
    """Finite-value mask, optionally excluding pixels at/above a saturation
    ceiling. Bit-depth-agnostic saturation detection isn't attempted beyond
    this explicit, caller-supplied threshold -- see lsf_stack.py's
    --saturation-level."""
    # Registered frames are zero-filled outside their footprint (dither
    # borders, rotation). Exact zeros are "no data", not signal -- Siril's
    # normalization/stacking also ignores them. Treating them as valid is
    # what produced dark, mismatched edges.
    mask = np.isfinite(frame) & (frame != 0)
    if saturation_level is not None:
        mask &= frame < saturation_level
    return mask


def _block_shape(h: int, w: int, target_blocks_per_axis: int = 48, min_block_size: int = 4) -> "tuple[int, int]":
    """Block height/width used by _block_medians (shared so the fit can work out how
    many pixels each block median averages)."""
    return (max(min_block_size, h // max(1, target_blocks_per_axis)),
            max(min_block_size, w // max(1, target_blocks_per_axis)))


def _block_medians(source: np.ndarray, reference: np.ndarray, valid: np.ndarray,
                    target_blocks_per_axis: int = 48, min_valid_frac: float = 0.25,
                    min_block_size: int = 4):
    """Tile (source, reference) into a target_blocks_per_axis x
    target_blocks_per_axis-ish grid and return the per-block median of each,
    for every block with enough valid pixels. Block *count* is bounded by
    target_blocks_per_axis regardless of image resolution (block *size*
    scales up instead), which is what keeps the O(n^2) Theil-Sen pass below
    tractable on any sensor size.

    `min_block_size` is a floor on block dimensions: for a real sensor
    frame (thousands of pixels per axis) target_blocks_per_axis alone
    already yields blocks comfortably above this floor, but for a small
    image -- target_blocks_per_axis=48 on a 64px-wide canonical grid would
    otherwise give floor(64/48)=1px blocks, i.e. no averaging at all,
    silently degrading back to the noise-dominated raw-pixel regime this
    function exists to avoid. The floor costs nothing on real data and
    keeps small grids (and small test fixtures) from quietly falling back
    to unaveraged pixels.

    Averaging within a block does two useful things at once: it suppresses
    per-pixel photometric noise (by roughly sqrt(block area)) relative to
    the smooth large-scale structure (gradients, vignetting) the fit is
    supposed to lock onto, and the median (rather than mean) naturally
    drops point-source contamination -- a star, hot pixel, or cosmic ray
    hit occupying a handful of pixels within an otherwise-background block
    simply doesn't move that block's median.
    """
    h, w = source.shape
    bs_y, bs_x = _block_shape(h, w, target_blocks_per_axis, min_block_size)
    bx, by = [], []
    for r0 in range(0, h, bs_y):
        for c0 in range(0, w, bs_x):
            vb = valid[r0:r0 + bs_y, c0:c0 + bs_x]
            if vb.sum() < min_valid_frac * vb.size:
                continue
            xb = source[r0:r0 + bs_y, c0:c0 + bs_x]
            yb = reference[r0:r0 + bs_y, c0:c0 + bs_x]
            bx.append(np.median(xb[vb]))
            by.append(np.median(yb[vb]))
    return np.asarray(bx, dtype=np.float64), np.asarray(by, dtype=np.float64)


def _theil_sen(bx: np.ndarray, by: np.ndarray) -> "tuple[float, float] | None":
    """Median-of-all-pairwise-slopes estimator: scale = median over every
    pair (i, j) of (by[j]-by[i])/(bx[j]-bx[i]), offset = median(by -
    scale*bx). Deliberately used INSTEAD OF an ordinary least-squares seed:
    OLS on the raw (even lightly cleaned) data is not robust to a handful
    of extreme-leverage points (an uncorrected satellite trail spanning
    many blocks at a wildly different value is exactly this), and once the
    very first fit is thrown off by them, residual-based iterative
    clipping built on top of that bad fit can fail to recover -- this was
    verified empirically to fail (scale collapsing toward 0) before this
    function replaced it. Theil-Sen's breakdown point (~29%) comfortably
    covers realistic contamination levels and has no such bootstrapping
    problem: every pairwise slope is an independent, self-contained
    estimate, so a bad seed can't propagate.

    Returns None if there are too few blocks (<8) to form a meaningful
    pairwise estimate.
    """
    n = bx.size
    if n < 8:
        return None
    ii, jj = np.triu_indices(n, k=1)
    dx = bx[jj] - bx[ii]
    span = np.percentile(bx, 90) - np.percentile(bx, 10)
    min_dx = 1e-9 if span == 0 else max(1e-12, 1e-6 * span)
    keep = np.abs(dx) > min_dx
    if keep.sum() < 8:
        return None
    slopes = (by[jj[keep]] - by[ii[keep]]) / dx[keep]
    scale = float(np.median(slopes))
    offset = float(np.median(by - scale * bx))
    return scale, offset


@dataclass
class NormFit:
    """Result of fit_normalization(): reference ~= source*scale + offset.

    method: "linear"      slope and offset both fitted (the normal case)
            "offset_only" the slope could not be determined (the reference has too little
                          large-scale structure, e.g. a flat field), so scale is 1 and only
                          an additive offset was matched
            "failed"      not enough overlap to fit anything; scale=1, offset=0 and the
                          caller should leave the frame unnormalized
    reason: human-readable cause for the two fallback methods (for logging)."""
    scale: float = 1.0
    offset: float = 0.0
    method: str = "failed"
    reason: str = ""



def _pixel_noise_sigma(frame: np.ndarray, valid: np.ndarray, max_rows: int = 1000) -> "float | None":
    """Robust per-pixel noise sigma from horizontal neighbour differences:
    1.4826 * MAD(d) / sqrt(2), d = x[i, j+1] - x[i, j]. Differencing removes the
    smooth sky/gradient, the MAD ignores stars and hot pixels, and only every
    k-th row is read so a huge frame costs a few MB. Used only for the ratio of
    two frames' noise (see _deming_slope), so a common bias such as correlated
    noise from resampling cancels. None if there is too little valid data."""
    step = max(1, frame.shape[0] // max_rows)
    f = np.asarray(frame[::step], dtype=np.float64)
    v = valid[::step]
    d = (f[:, 1:] - f[:, :-1])[v[:, 1:] & v[:, :-1]]
    if d.size < 1000:
        return None
    med = float(np.median(d))
    sig = 1.4826 * float(np.median(np.abs(d - med))) / np.sqrt(2.0)
    return sig if (np.isfinite(sig) and sig > 0) else None


def _deming_slope(x: np.ndarray, y: np.ndarray, delta: float) -> "tuple[float, float] | None":
    """Deming (errors-in-variables) regression of y on x, y = a*x + b, where BOTH
    x and y carry noise and delta = var(error in y) / var(error in x), the two
    error variances in their own units. delta = 1 is orthogonal regression;
    delta -> infinity recovers ordinary least squares of y on x. Unlike OLS it
    is not biased toward zero by noise in x."""
    n = x.size
    if n < 3 or not np.isfinite(delta) or delta <= 0:
        return None
    mx, my = float(x.mean()), float(y.mean())
    sxx = float(np.sum((x - mx) ** 2) / (n - 1))
    syy = float(np.sum((y - my) ** 2) / (n - 1))
    sxy = float(np.sum((x - mx) * (y - my)) / (n - 1))
    if sxy == 0:
        return None
    t = syy - delta * sxx
    a = (t + np.sqrt(t * t + 4.0 * delta * sxy * sxy)) / (2.0 * sxy)
    return float(a), float(my - a * mx)


def _slope_determined(x: np.ndarray, resid: np.ndarray, a: float, max_rel_err: float) -> "tuple[bool, str]":
    """Is the fitted slope statistically determined? The standard error of a
    least-squares slope is residual scatter / (sqrt(n) * spread of X); with too
    little structure in the reference the slope is noise. Requires a positive
    slope and standard error <= max_rel_err * slope. `x` and `resid` are the
    kept blocks' reference values and residuals about the line."""
    n_k = int(x.size)
    sxx = float(np.sum((x - x.mean()) ** 2)) if n_k else 0.0
    if not (np.isfinite(a) and a > 1e-6) or sxx <= 0 or n_k < 8:
        return False, "slope not positive (reference has too little large-scale structure to determine it)"
    s_res = float(np.sqrt(np.sum((resid - resid.mean()) ** 2) / max(n_k - 2, 1)))
    if (s_res / np.sqrt(sxx)) / a > max_rel_err:
        return False, "slope not determined (reference has too little large-scale structure)"
    return True, ""


def fit_normalization(source: np.ndarray, reference: np.ndarray,
                      valid_mask: "np.ndarray | None" = None,
                      reference_is_clean: bool = True,
                      reference_to_raw: "tuple[float, float] | None" = None,
                      target_blocks_per_axis: int = 48,
                      min_valid_frac: float = 0.25,
                      polish_sigma_clip: float = 5.0,
                      max_slope_rel_err: float = 0.1,
                      max_slope_disagreement: float = 2.0,
                      max_deming_gain: float = 2.0,
                      min_reference_reliability: float = 0.5,
                      trim_fraction: float = 0.05) -> NormFit:
    """
    Robust linear fit of `source` onto `reference` at matching pixel
    positions: reference ~= source*scale + offset. Implements
    live_stack.md Section 8.3's ingestion-time normalization as a
    block-median-then-Theil-Sen regression (see _block_medians /
    _theil_sen above for why each stage is there) followed by one light
    sigma-clipped polish pass on the block medians, for the edge case of
    contamination severe enough to still bias the initial Theil-Sen
    estimate a little.

    Direction of the regression (regression dilution): the line is fitted as
    source = a*reference + b, with the REFERENCE as the explanatory variable X
    and the noisy incoming frame as Y, and then inverted: scale = 1/a,
    offset = -b/a. Regressing the other way round (reference on source) puts
    the noisy variable on X, which biases the slope toward zero by
    var(structure) / (var(structure) + var(noise)); with a small gradient
    range the scale collapsed (measured: -34% at a gradient of 0.02 and -90%
    at 0.005, with 0.01 per-pixel noise on 10x10 px blocks). Noise in Y does
    not bias the slope.

    reference_is_clean: True when the reference is the accumulator's derived
    mean, which is much cleaner than one frame, so ordinary least squares with
    X = reference is (nearly) unbiased. False on first initialization, when the
    reference is a single batch frame: both variables are then noisy, and OLS
    leaves a residual bias (measured +4% in scale at a gradient of 0.02), so the
    final slope is a Deming errors-in-variables fit instead, with the error
    variance ratio taken from each frame's own measured noise
    (_pixel_noise_sigma). The Theil-Sen stage and the outlier clip are
    unchanged; Deming only replaces the final least-squares step.

    reference_to_raw: (gain, bias) such that raw_reference = gain*reference + bias.
    Used when the reference plane is stored output-normalized (Stack plane min-max
    mapped to [0, 1]): the block medians are converted to raw units before fitting
    (an affine map commutes with the median), so scale and offset come out directly
    in raw accumulator units and no after-the-fact conversion is needed. This is
    also what makes the additive fallback below correct: "scale = 1" must mean 1 in
    raw units.

    Fallbacks. Returns a NormFit:
      - "failed" (scale=1, offset=0): fewer than 8 valid blocks, i.e. not enough
        overlap to fit anything. The caller should leave the frame unnormalized.
      - "offset_only": there is overlap but the slope is not statistically
        determined, i.e. the reference has too little large-scale structure for
        the slope's standard error to be below `max_slope_rel_err` of the slope
        (a nearly flat field), or the slope is not positive. A number for the
        scale would be noise, so scale = 1 and only an additive offset is matched
        (median of reference - source over the blocks, clipped at 5 MAD), which is
        the one quantity a flat field does determine.
      - "linear": the normal case.

    A slope is also rejected as not determined when the robust seed and the refined
    slope disagree by more than `max_slope_disagreement` (see the stability check), and,
    on the first-init (Deming) path, when the errors-in-variables correction exceeds
    `max_deming_gain` times the least-squares slope, and, on the clean-reference path,
    when the reference's own block-median noise leaves a reliability below
    `min_reference_reliability` (the reference is then not clean enough for plain least
    squares). Healthy fits measure a reliability near 1. Finally the slope must be
    determined on the central blocks alone (the `trim_fraction` most extreme reference
    blocks dropped on each side), so that it is supported by the bulk of the field and not
    by a few bright-star blocks on an otherwise structureless sky.

    Known limitation: like Siril's -norm=addscale or PixInsight's LinearFit used
    between two images, a slope needs *some* shared large-scale structure
    (gradient, vignetting, nebulosity) across a reasonable fraction of the
    field. A field that is almost perfectly flat except for a sparse handful of
    stars cannot give a trustworthy scale; this is a property of the technique
    itself. In that case the fallback matches the sky level only, and any real
    transparency change between the two frames is left uncorrected.
    """
    mask = np.isfinite(source) & np.isfinite(reference)
    if valid_mask is not None:
        mask &= valid_mask

    bsrc, bref = _block_medians(source, reference, mask, target_blocks_per_axis, min_valid_frac)
    if reference_to_raw is not None:
        bref = bref * float(reference_to_raw[0]) + float(reference_to_raw[1])
    if bsrc.size < 8:
        return NormFit(1.0, 0.0, "failed", "fewer than 8 valid overlapping blocks")

    reason = ""
    linear = None
    # X = reference, Y = noisy source  (see docstring: regression dilution)
    fit = _theil_sen(bref, bsrc)
    if fit is None:
        reason = "reference has no usable large-scale structure"
    else:
        a, b = fit                                   # source = a*reference + b
        a_seed = a

        # Light polish: drop any block that's still a clear outlier under the
        # Theil-Sen fit, then refit on the rest.
        resid = bsrc - (bref * a + b)
        med = float(np.median(resid))
        mad = float(np.median(np.abs(resid - med)) * 1.4826 + 1e-12)
        keep = np.abs(resid - med) <= polish_sigma_clip * mad
        if keep.sum() < 8:
            keep = np.ones(bref.size, dtype=bool)

        xk = bref[keep]
        a_ls = a
        if reference_is_clean:
            # Unchanged accumulator path: Theil-Sen line, or least squares after the clip.
            if not np.all(keep):             # (keep always holds >= 8 blocks here)
                A = np.vstack([bref[keep], np.ones(keep.sum())]).T
                (a, b), *_ = np.linalg.lstsq(A, bsrc[keep], rcond=None)
                resid = bsrc - (bref * a + b)
                a_ls = a
        else:
            # Both axes noisy (first init). Decide FIRST whether any slope is real, using
            # the ordinary least-squares slope of source on reference: its standard error
            # formula is valid, and for a field whose blocks form a structureless cloud it
            # correctly fails. (The Deming slope is not safe to test this way: on a
            # structureless cloud it returns an arbitrary large slope that looks "precise".)
            A = np.vstack([xk, np.ones(xk.size)]).T
            (a_o, b_o), *_ = np.linalg.lstsq(A, bsrc[keep], rcond=None)
            a_ls = a_o
            resid_o = bsrc - (bref * a_o + b_o)
            ok_o, why_o = _slope_determined(xk, resid_o[keep], a_o, max_slope_rel_err)
            if not ok_o:
                a, reason = float("nan"), why_o       # force the fallback below
            else:
                a, b, resid = a_o, b_o, resid_o       # least squares unless Deming applies
                sig_s = _pixel_noise_sigma(source, mask)
                sig_r = _pixel_noise_sigma(reference, mask)
                if sig_s is not None and sig_r is not None:
                    if reference_to_raw is not None:
                        sig_r *= abs(float(reference_to_raw[0]))
                    dm = _deming_slope(xk, bsrc[keep], (sig_s / sig_r) ** 2)
                    if dm is not None and dm[0] > 0:
                        # Deming/OLS slope ~ 1/reliability, where reliability = the share of the
                        # reference's block-to-block variance that is real structure rather than
                        # noise. A large correction means noise dominates the structure and the
                        # corrected slope, though plausible-looking, is poorly determined.
                        if dm[0] / a_o > max_deming_gain:
                            a = float("nan")
                            reason = ("noise dominates the structure in the reference (errors-in-variables "
                                      f"correction {dm[0] / a_o:.1f}x is too large to trust)")
                        else:
                            a, b = dm
                            resid = bsrc - (bref * a + b)

        if np.isfinite(a):
            ok_s, why = _slope_determined(xk, resid[keep], a, max_slope_rel_err)
            # Stability: the robust Theil-Sen seed and the least-squares slope after the
            # outlier clip must agree. On a field that is mostly structureless with a few
            # bright features (stars on a flat sky) the seed is just noise and the clip,
            # anchored on it, throws away the informative blocks; the surviving line can
            # still pass the standard-error test and be confidently wrong. Measured on
            # healthy fits the two never differ by more than 1.2x; in that failure ~6x.
            if ok_s:
                # Bulk support: the slope must also be determined on the central blocks alone
                # (dropping the `trim_fraction` most extreme reference blocks on each side).
                # On a nearly flat field the informative blocks are a few percent of the total
                # (bright stars); the bulk is a structureless cloud, so the surviving line
                # is not a gradient/vignetting fit at all, just noise plus a few leverage points.
                lo_t, hi_t = np.percentile(xk, [100 * trim_fraction, 100 * (1 - trim_fraction)])
                inner = (xk >= lo_t) & (xk <= hi_t)
                if inner.sum() >= 8:
                    xi, yi = xk[inner], bsrc[keep][inner]
                    Ai = np.vstack([xi, np.ones(xi.size)]).T
                    (a_i, b_i), *_ = np.linalg.lstsq(Ai, yi, rcond=None)
                    ok_i, _why_i = _slope_determined(xi, yi - (xi * a_i + b_i), a_i, max_slope_rel_err)
                    if not ok_i:
                        ok_s = False
                        why = ("slope not supported by the bulk of the field (it rests on a few "
                               "extreme blocks such as bright stars; the rest is structureless)")
            if ok_s and reference_is_clean:
                # Reliability of the reference: the share of its block-to-block variance that
                # is real structure rather than noise in its block medians. Clean-path least
                # squares attenuates the slope by this factor, so it is only (nearly) unbiased
                # when the reference is much cleaner than the structure it is used for. A
                # young accumulator over a nearly flat field is not (cloud of blocks plus a
                # few faint stars: the fitted slope came out ~0.3 where the truth was ~1).
                sig_r = _pixel_noise_sigma(reference, mask)
                if sig_r is not None and xk.size > 2:
                    if reference_to_raw is not None:
                        sig_r *= abs(float(reference_to_raw[0]))
                    bs_y, bs_x = _block_shape(reference.shape[0], reference.shape[1], target_blocks_per_axis)
                    var_noise = (np.pi / 2.0) * sig_r * sig_r / float(bs_y * bs_x)   # variance of a block median
                    var_x = float(np.var(xk, ddof=1))
                    rel = max(var_x - var_noise, 0.0) / var_x if var_x > 0 else 0.0
                    if rel < min_reference_reliability:
                        ok_s = False
                        why = (f"reference too noisy for its structure (reliability {rel:.2f} < "
                               f"{min_reference_reliability:g}): least squares would shrink the slope")
            if ok_s:
                lo_s, hi_s = min(a_seed, a_ls), max(a_seed, a_ls)
                if not (lo_s > 0 and hi_s / lo_s <= max_slope_disagreement):
                    ok_s = False
                    why = ("slope unstable (robust and refined slopes disagree: the field is mostly "
                           "structureless with a few bright features)")
            if ok_s:
                linear = NormFit(float(1.0 / a), float(-b / a), "linear")     # invert: ref = src/a - b/a
            else:
                reason = why
    if linear is not None:
        return linear

    d = bref - bsrc
    med = float(np.median(d))
    mad = float(np.median(np.abs(d - med)) * 1.4826 + 1e-12)
    ok = np.abs(d - med) <= 5.0 * mad
    return NormFit(1.0, float(np.median(d[ok])) if ok.any() else med, "offset_only", reason)


NORM_METHODS = ("linearfit", "addscale")
NORM_ESTIMATORS = ("lite", "siril")


def frame_estimators(plane: np.ndarray, valid: np.ndarray, estimator: str = "lite",
                     stats=None) -> "tuple[float, float] | None":
    """(location, scale) of one colour plane of one frame, for the 'addscale' method.

    estimator "lite":  median and 1.5 * MAD over the valid (finite, non-zero) pixels. This is
                       what Siril's -lite_norm uses (1.5*MAD approximates sqrt(bwmv)).
    estimator "siril": the IKSS location / scale Siril measured on the frame (`stats` is the
                       sirilpy statistics object: .location, .scale). They are measured on the
                       frame as Siril held it when the metadata was gathered, i.e. BEFORE any
                       alignment resampling.
    Returns None when the estimate is unavailable or not finite."""
    if estimator == "siril":
        if stats is None:
            return None
        try:
            loc, sca = float(stats.location), float(stats.scale)
        except (AttributeError, TypeError, ValueError):
            return None
        return (loc, sca) if (np.isfinite(loc) and np.isfinite(sca)) else None
    if estimator != "lite":
        raise ValueError(f"unknown normalization estimator '{estimator}' (supported: {NORM_ESTIMATORS})")
    v = np.asarray(plane)[valid]
    if v.size < 16:
        return None
    v = v.astype(np.float64, copy=False)
    loc = float(np.median(v))
    mad = float(np.median(np.abs(v - loc)))
    return loc, 1.5 * mad


def addscale_fit(loc: float, scale: float, loc_ref: float, scale_ref: float) -> NormFit:
    """Siril's -norm=addscale for one frame and channel, in LSF's convention (x*scale + offset).

    Siril (normalization.c, compute_factors_from_estimators + median_and_mean.c):
        s = scale_ref / scale  (1 if scale == 0);  poffset = s*loc - loc_ref;  pixel' = pixel*s - poffset
    which is  pixel' = s*(pixel - loc) + loc_ref,  i.e.  offset = loc_ref - s*loc  here."""
    if not all(np.isfinite(v) for v in (loc, scale, loc_ref, scale_ref)):
        return NormFit(1.0, 0.0, "failed", "non-finite location/scale estimate")
    s = 1.0 if scale == 0 else float(scale_ref / scale)
    return NormFit(s, float(loc_ref - s * loc), "addscale")


def apply_linear_normalization(frame: np.ndarray, scale: float, offset: float) -> np.ndarray:
    """reference ~= source*scale + offset, applied to the frame's full
    extent (not just whatever overlap region the coefficients were fit
    from) -- so pixels the accumulator hasn't seen before (e.g. a dither
    extending the canonical grid's footprint, Section 7) still get merged
    in using the same photometric normalization derived from the
    overlap that does exist."""
    return frame * scale + offset


# --------------------------------------------------------------------------
# Batch statistics: iterative, per-pixel rejection across a stack of frames
# held fully in memory, producing a weighted mean / M2 / weight / count per
# pixel -- i.e. exactly what a single ordinary Siril "stack ... rej ..."
# call would produce, computed here so it can be merged into (or become) an
# LSF accumulator.
# --------------------------------------------------------------------------

SUPPORTED_REJECTION = ("none", "sigma", "mad", "winsorized", "median_sigma", 
                       "percentile", "linearfit", "gesdt")
SUPPORTED_WEIGHTING = ("none", "noise", "nbstars", "wfwhm")

# Per-method (sigma_low, sigma_high) defaults, following Siril's own defaults.
# percentile thresholds are FRACTIONS of the median (0.2 = 20%), so a generic
# 3.0 default is meaningless there. gesdt ignores these: it uses gesdt_max_frac
# (max outlier fraction) and gesdt_alpha (significance) instead.
REJECTION_DEFAULTS = {
    "none": (3.0, 3.0), "sigma": (3.0, 3.0), "mad": (3.0, 3.0),
    "winsorized": (3.0, 3.0), "median_sigma": (3.0, 3.0),
    "linearfit": (5.0, 5.0), "percentile": (0.2, 0.1), "gesdt": (3.0, 3.0),
}


def resolve_rejection_params(method, sigma_low, sigma_high):
    dl, dh = REJECTION_DEFAULTS[method]
    return (dl if sigma_low is None else float(sigma_low),
            dh if sigma_high is None else float(sigma_high))

# Rejection algorithms: All operate on the current batch held in memory.
# They do NOT operate across accumulated history (that would violate the
# bounded-accumulator design). All methods are modelled on Siril's, not
# reproduced line for line; the docstring of compute_batch_stats_numpy lists
# the known differences (rejection order at the 4-sample limit, passes vs
# convergence, sigma estimator). They suit the Pass 2 rejection step within
# lsf_stack.py where the full batch is resident.
# True iterative sigma clipping across the whole history remains unimplemented
# as it is fundamentally incompatible with a bounded accumulator (Section 5.2).


# --------------------------------------------------------------------------
# Batch statistics: iterative, per-pixel rejection across a stack of frames
# --------------------------------------------------------------------------

# Minimum number of samples a pixel keeps after iterative rejection. Siril stops
# rejecting once only 4 samples remain ("N - r <= 4" is checked before each
# further rejection), so up to (k - 4) of a pixel's k active samples can go.
MIN_KEEP = 4


def _limit_bad_numpy(bad, k, dev, sigma, sigma_low, sigma_high, min_keep=MIN_KEEP):
    """Cap the rejections per pixel at max(k - min_keep, 0), where k is the
    number of samples currently active at that pixel, keeping the most
    extreme candidates and dropping the rest of the candidates for this pass.

    The previous behaviour cancelled EVERY rejection at a pixel as soon as the
    full candidate set would have left too few samples, so a pixel with one
    obvious outlier plus one marginal candidate kept both. Ranking by
    extremity instead stops exactly at the limit.

    Extremity = |dev| in units of its (asymmetric) clipping threshold,
    sigma_low*sigma below / sigma_high*sigma above. It is computed only for
    pixels that are over budget, so the common case costs one comparison per
    pixel.

    bad, dev: (N, ...) ; k, sigma: (...) active-sample count and sigma per pixel.
    """
    budget = np.maximum(k - min_keep, 0)
    over = bad.sum(axis=0) > budget
    if not over.any():
        return bad
    n = bad.shape[0]
    pix = np.nonzero(over)
    cols = (slice(None),) + pix
    sub_bad = bad[cols]                                    # (N, M)
    d = dev[cols]
    with np.errstate(divide="ignore", invalid="ignore"):
        score = np.abs(d) / (np.where(d < 0, sigma_low, sigma_high) * sigma[pix][None, :])
    sub_s = np.where(sub_bad, score, -np.inf)
    order = np.argsort(-sub_s, axis=0, kind="stable")      # most extreme first
    rank = np.empty_like(order)
    np.put_along_axis(rank, order, np.arange(n)[:, None], axis=0)
    out = bad.copy()
    out[cols] = sub_bad & (rank < budget[over][None, :])
    return out


def _active_sample_sd(stack, active):
    """Sample standard deviation (N-1 denominator), UNWEIGHTED, about the mean
    of the active samples, per pixel: the estimator Siril's sigma clipping uses
    (siril_stats_float_sd: two-pass, N-1, no weights). Pixels with fewer than 2
    active samples get 0. `stack` is already zero where inactive."""
    k = active.sum(axis=0)
    m = (stack * active).sum(axis=0) / np.maximum(k, 1)
    d = np.where(active, stack - m, 0.0)
    return np.sqrt((d * d).sum(axis=0) / np.maximum(k - 1, 1))


def _linearfit_reject_numpy(stack, active, sigma_low, sigma_high):
    """Vectorized linear-fit clipping (Siril-style), all pixels at once.

    Per pixel: sort the active samples, fit a line value ~ rank, take the mean
    absolute residual as sigma, reject samples beyond sigma_low/sigma_high of
    the line, and repeat until a pass rejects nothing. Each pixel never drops
    below MIN_KEEP (4) samples -- when that limit binds, the most extreme
    candidates are rejected first -- and pixels with <= 4 active samples are left alone, so
    empty border pixels (count 0) no longer affect convergence of the others.
    Same algorithm as the torch backend.
    """
    n = stack.shape[0]
    bshape = (n,) + (1,) * (stack.ndim - 1)
    idx = np.arange(n, dtype=np.float64).reshape(bshape)
    active = active.copy()
    for _ in range(max(n, 1)):
        k = active.sum(axis=0)
        if not np.any(k > 4):
            break
        xs_in = np.where(active, stack, np.inf)
        order = np.argsort(xs_in, axis=0, kind="stable")
        xs = np.take_along_axis(xs_in, order, axis=0)
        valid = idx < k[None]
        ys = np.where(valid, xs, 0.0).astype(np.float64)
        kf = np.maximum(k, 1).astype(np.float64)
        m_x = (kf - 1.0) / 2.0
        m_y = ys.sum(axis=0) / kf
        di = idx - m_x[None]
        dy = np.where(valid, ys - m_y[None], 0.0)
        var = np.where(valid, di * di, 0.0).sum(axis=0) / kf
        slope = ((di * dy).sum(axis=0) / kf) / np.maximum(var, 1e-300)
        icpt = m_y - m_x * slope
        expected = icpt[None] + slope[None] * idx
        resid = ys - expected
        sig = np.where(valid, np.abs(resid), 0.0).sum(axis=0) / kf
        cand = valid & ((-resid > sig[None] * sigma_low) | (resid > sig[None] * sigma_high)) \
            & (k > 4)[None]
        # Positions here are in sorted order; keep the MIN_KEEP samples closest to the
        # line and, when over budget, reject the most extreme candidates first.
        allowed = _limit_bad_numpy(cand, k, resid, sig, sigma_low, sigma_high)
        if not allowed.any():
            break
        rej = np.zeros(allowed.shape, dtype=bool)
        np.put_along_axis(rej, order, allowed, axis=0)      # back to original frame order
        active &= ~rej
    return active


def _gesdt_block(x, act, table, max_frac):
    """Generalized ESD (Rosner) rejection for a block of pixels: x, act are (n, P).
    Returns the new active mask (n, P).

    Per pixel: sort the active samples; for up to floor(k * max_frac) rounds, take the
    window's extreme sample (the one farthest from the window mean, high end on a tie
    of the two ends going to the low end), compute G = deviation / sample SD (N-1) and
    compare with the critical value for the current window size; then drop that sample
    from the window. Rosner's rule: everything up to the LAST significant candidate is
    an outlier. The window is always a contiguous range [lo, hi] of the sorted samples,
    so the window sum and sum of squares are updated by subtracting the removed sample
    instead of being recomputed (O(1) per round instead of O(n)). Values are shifted
    by the pixel's mean first so the running sum of squares does not cancel badly.
    """
    n, P = x.shape
    k = act.sum(axis=0)
    kf = k.astype(np.float64)
    mo = np.maximum(np.floor(kf * max_frac), 1.0)
    mo = np.where(mo >= kf, kf - 1.0, mo)
    mo = np.where(k > 1, mo, 0.0).astype(np.int64)
    T = int(mo.max()) if P else 0
    if T <= 0:
        return act

    xin = np.where(act, x.astype(np.float64), np.inf)
    order = np.argsort(xin, axis=0, kind="stable")
    ys = np.take_along_axis(xin, order, axis=0)
    valid = np.arange(n)[:, None] < k[None, :]
    ctr = np.where(valid, ys, 0.0).sum(axis=0) / np.maximum(k, 1)
    d = np.where(valid, ys - ctr[None, :], 0.0)           # shifted sorted samples
    del xin, ys, valid
    S = d.sum(axis=0)
    Q = (d * d).sum(axis=0)
    lo = np.zeros(P, dtype=np.int64)
    hi = np.maximum(k - 1, 0).astype(np.int64)
    cols = np.arange(P)
    pos = np.empty((T, P), dtype=np.int32)                # sorted position tested at round t
    last = np.full(P, -1, dtype=np.int64)                 # last significant round per pixel
    std_floor = 1e-9 * np.abs(ctr) + 1e-300               # SD below this is "no scatter"
    for t in range(T):
        vt = mo > t
        size = hi - lo + 1
        szf = np.maximum(size, 1).astype(np.float64)
        mean_s = S / szf
        var = np.maximum(Q - szf * mean_s * mean_s, 0.0) / np.maximum(size - 1, 1)
        std = np.sqrt(var)
        d_lo = d[lo, cols]
        d_hi = d[hi, cols]
        dlow = mean_s - d_lo
        dhigh = d_hi - mean_s
        use_high = dhigh > dlow
        dev = np.where(use_high, dhigh, dlow)
        good = std > std_floor
        G = np.where(good, dev / np.where(good, std, 1.0), 0.0)
        sig = vt & (G > table[np.clip(size, 0, n)])
        last = np.where(sig, t, last)
        pos[t] = np.where(use_high, hi, lo)
        rem = np.where(use_high, d_hi, d_lo)
        S = S - np.where(vt, rem, 0.0)
        Q = Q - np.where(vt, rem * rem, 0.0)
        lo = lo + (vt & ~use_high)
        hi = hi - (vt & use_high)
    if not (last >= 0).any():
        return act
    rej_sorted = np.zeros((n, P), dtype=bool)
    for t in range(T):
        m = last >= t                                      # t <= last significant round
        if not m.any():
            break
        rej_sorted[pos[t][m], cols[m]] = True
    rej = np.zeros((n, P), dtype=bool)
    np.put_along_axis(rej, order, rej_sorted, axis=0)      # back to original frame order
    return act & ~rej


def _gesdt_reject_numpy(stack, active, gesdt_alpha, gesdt_max_frac, block_samples=2_000_000):
    """Vectorized GESDT rejection over a strip: stack/active are (n, ...) arrays.
    Pixels are processed in column blocks of about `block_samples` samples so the
    temporaries stay small whatever the strip size. Returns a new active mask."""
    if not SCIPY_AVAILABLE:
        raise RuntimeError("GESDT rejection requires scipy. Install with: pip install scipy")
    n = stack.shape[0]
    shape = stack.shape
    P = int(np.prod(shape[1:], dtype=np.int64))
    x2 = np.asarray(stack).reshape(n, P)
    a2 = np.asarray(active).reshape(n, P)
    if not (a2.sum(axis=0) > 1).any():
        return active
    table = _esd_crit_table(n, gesdt_alpha)
    out = a2.copy()
    chunk = max(1024, block_samples // max(n, 1))
    for c0 in range(0, P, chunk):
        sl = slice(c0, min(P, c0 + chunk))
        out[:, sl] = _gesdt_block(x2[:, sl], a2[:, sl], table, gesdt_max_frac)
    return out.reshape(shape)


def compute_batch_stats_numpy(
    stack: np.ndarray,          # shape (N, H, W) or (N, H, W, C), already normalized
    frame_weights: np.ndarray,  # shape (N,)
    rejection: str = "sigma",
    sigma_low: float = 3.0,
    sigma_high: float = 3.0,
    max_iters: int = 3,
    gesdt_alpha: float = 0.05,
    gesdt_max_frac: float = 0.3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute per-pixel weighted mean, M2 (Welford sum of squared weighted
    deviations), total weight, sum of squared weights and surviving-frame count
    for one batch.

    Returns (mean, m2, weight, w2, count), each shaped like stack.shape[1:].
    `w2` is the sum over the accepted samples of (frame weight)^2 -- the LSF_W2
    plane of the FITS profile (Section 3). Together with `weight` it gives the
    effective sample size n_eff = weight^2 / w2 (Section 10), which, unlike
    `count`, is not affected by the uint16 saturation of the count plane.
    
    Modelled on Siril's rejection (rejection_float.c), not a line-for-line
    reproduction. Known differences from Siril are noted at each method:
    - iterative methods (sigma, mad, winsorized) stop rejecting at a pixel
      once MIN_KEEP (4) samples remain, as Siril does, but when the limit
      binds they reject the MOST EXTREME candidates first; Siril rejects in
      frame order;
    - `max_iters` counts rejection passes (Siril iterates to convergence);
    - rejection ignores frame weights (Siril's sigma estimator is unweighted);
      weights apply only to the final mean;
    - "median_sigma" uses the same test as "sigma" (Siril documents it as the same
      algorithm with rejected samples replaced by the median instead of excluded),
      without the zero-median guard and without the 4-sample limit (replacement
      never reduces the sample count);
    """
    if rejection not in SUPPORTED_REJECTION:
        raise ValueError(f"unsupported rejection method '{rejection}' (supported: {SUPPORTED_REJECTION})")

    n = stack.shape[0]
    w = frame_weights.reshape((n,) + (1,) * (stack.ndim - 1)).astype(np.float64)
    # NaN marks "this frame has no data at this pixel" (outside footprint).
    # Such samples start inactive so they never get weight, count, or
    # influence rejection statistics.
    active = np.isfinite(stack)
    stack = np.where(active, stack, 0.0)

    # Phase 1: Pre-compute pixel-wise eligibility masks (matching Siril's early exits)
    # Siril rejection_float.c line 157: "if (kept <= 1) return kept"
    active_count = active.sum(axis=0)
    can_attempt_rejection = active_count > 1
    
    # `iters` is the number of rejection PASSES (none for "none" or fewer than 3 frames).
    # The final weighted statistics are recomputed after the loop, so the last pass needs
    # no extra mean update.
    iters = 0 if (rejection == "none" or n < 3) else max_iters

    for it in range(iters):
        if rejection == "mad":
            # Robust estimator throughout: median location, MAD-based scale.
            # Siril rejection_float.c lines 168-170: zero median check
            # Siril rejection_float.c lines 192-225: MAD rejection loop
            
            # Suppress NaN warnings for pixels outside footprint
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'All-NaN slice')
                med = np.nanmedian(np.where(active, stack, np.nan), axis=0)
            
            # Siril lines 168-170: "if (median == 0.0) return 0"
            median_is_valid = (med != 0.0) & np.isfinite(med)
            can_reject_this_iter = can_attempt_rejection & median_is_valid
            
            if not np.any(can_reject_this_iter):
                break  # No pixels eligible for rejection
            
            abs_dev = np.abs(stack - med)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'All-NaN slice')
                sigma = np.nanmedian(np.where(active, abs_dev, np.nan), axis=0) * 1.4826 + 1e-12
            
            dev = stack - med
            lo = -sigma_low * sigma
            hi = sigma_high * sigma
            newly_bad = active & ((dev < lo) | (dev > hi))
            
            # Only reject on eligible pixels
            newly_bad = newly_bad & can_reject_this_iter
            
            # Never leave a pixel with fewer than MIN_KEEP samples; if that limit
            # binds, reject the most extreme candidates first (rest are kept).
            newly_bad = _limit_bad_numpy(newly_bad, active.sum(axis=0), dev, sigma,
                                         sigma_low, sigma_high)
            
            # Exclude rejected samples (a pixel can no longer be emptied: at most
            # k - MIN_KEEP of its k active samples are rejected).
            still_active = active & (~newly_bad)
            if not np.array_equal(still_active, active):
                active = still_active
            else:
                break
                
        elif rejection == "sigma":
            # Siril rejection_float.c lines 168-170: zero median check
            # Siril rejection_float.c lines 191-225: sigma rejection loop
            
            # Compute median for zero-check (Siril lines 168-170)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'All-NaN slice')
                med_check = np.nanmedian(np.where(active, stack, np.nan), axis=0)
            
            # Siril lines 168-170: "if (median == 0.0) return 0"
            median_is_valid = (med_check != 0.0) & np.isfinite(med_check)
            can_reject_this_iter = can_attempt_rejection & median_is_valid
            
            if not np.any(can_reject_this_iter):
                break  # No pixels eligible for rejection
            
            # Same estimator and centre as Siril's sigma clipping: the clip test is
            # taken about the MEDIAN of the active samples, in units of the
            # unweighted sample standard deviation (N-1) about their MEAN.
            # Frame weights do not enter the rejection decision, only the final mean.
            # (Established from Siril's documentation and a third-party port of
            # siril_stats_float_sd; not checked against rejection_float.c itself.)
            dev = stack - med_check
            sigma = _active_sample_sd(stack, active) + 1e-12
            
            lo = -sigma_low * sigma
            hi = sigma_high * sigma
            newly_bad = active & ((dev < lo) | (dev > hi))
            
            # Only reject on eligible pixels
            newly_bad = newly_bad & can_reject_this_iter
            
            # Never leave a pixel with fewer than MIN_KEEP samples; if that limit
            # binds, reject the most extreme candidates first (rest are kept).
            newly_bad = _limit_bad_numpy(newly_bad, active.sum(axis=0), dev, sigma,
                                         sigma_low, sigma_high)
            
            # Exclude rejected samples (a pixel can no longer be emptied).
            still_active = active & (~newly_bad)
            if not np.array_equal(still_active, active):
                active = still_active
            else:
                break
                
        elif rejection == "winsorized":
            # Winsorized sigma clipping, modelled on Siril's algorithm (1.5 sigma
            # winsorization, 1.134 correction, converge to 0.0005). Not verified
            # line-for-line against Siril, so treat as "close to", not "identical".
            # Siril rejection_float.c lines 240-276
            # Note: WINSORIZED does NOT have the zero-median check
            # (Siril lines 168-170 only apply to PERCENTILE/SIGMA/MAD)
            
            # Only attempt rejection where we have >1 active frames
            if not np.any(can_attempt_rejection):
                break
            
            # Suppress warnings for pixels with all-NaN slices (outside footprint)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'All-NaN slice')
                medianf = np.nanmedian(np.where(active, stack, np.nan), axis=0)
            
            # Handle NaN medians (pixels with no active frames)
            median_is_finite = np.isfinite(medianf)
            can_reject_winsor = can_attempt_rejection & median_is_finite
            
            if not np.any(can_reject_winsor):
                break
            
            # Weighted mean of the currently active samples (only this method needs it).
            wa = np.where(active, w, 0.0)
            total_w = wa.sum(axis=0)
            total_w_safe = np.where(total_w > 0, total_w, 1.0)
            mean = (wa * stack).sum(axis=0) / total_w_safe

            # Stage 1: Iteratively compute Winsorized sigma
            # Siril lines 243-254: initial sigma and inner convergence loop
            dev_ws = stack - mean
            var_ws = (wa * dev_ws * dev_ws).sum(axis=0) / total_w_safe
            sigma_ws = np.sqrt(np.maximum(var_ws, 0.0)) + 1e-12
            
            # Create working copy for Winsorization (Siril line 245)
            w_stack = stack.copy()
            
            # Inner convergence loop (Siril's do-while at lines 246-254)
            max_inner_iters = 20
            for inner_it in range(max_inner_iters):
                sigma_ws_prev = sigma_ws.copy()
                
                # Winsorize at ±1.5σ around median (Siril lines 247-250)
                m0 = medianf - 1.5 * sigma_ws
                m1 = medianf + 1.5 * sigma_ws
                
                # Clip w_stack values to [m0, m1]
                w_stack = np.clip(w_stack, m0, m1)
                
                # Compute weighted std dev of clipped stack with 1.134 correction (Siril line 253)
                dev_winsor = w_stack - mean
                var_winsor = (wa * dev_winsor * dev_winsor).sum(axis=0) / total_w_safe
                sigma_ws = 1.134 * np.sqrt(np.maximum(var_winsor, 0.0)) + 1e-12
                
                # Check convergence (Siril line 254): |sigma - sigma_prev| / sigma_prev < 0.0005
                rel_change = np.abs(sigma_ws - sigma_ws_prev) / (sigma_ws_prev + 1e-12)
                if np.all(rel_change < 0.0005):
                    break
            
            # Stage 2: Use converged Winsorized sigma for actual rejection
            # Siril lines 255-263: rejection with user's thresholds
            dev = stack - medianf
            lo_winsor = -sigma_low * sigma_ws
            hi_winsor = sigma_high * sigma_ws
            newly_bad = active & ((dev < lo_winsor) | (dev > hi_winsor))
            
            # Only reject on eligible pixels
            newly_bad = newly_bad & can_reject_winsor
            
            # Never leave a pixel with fewer than MIN_KEEP samples; if that limit
            # binds, reject the most extreme candidates first (rest are kept).
            newly_bad = _limit_bad_numpy(newly_bad, active.sum(axis=0), dev, sigma_ws,
                                         sigma_low, sigma_high)
            
            # Actually exclude rejected samples (Siril lines 266-274)
            still_active = active & (~newly_bad)
            if not np.array_equal(still_active, active):
                active = still_active
            else:
                break

        elif rejection == "median_sigma":
            # Siril documents Median Sigma Clipping as Sigma Clipping in which each
            # rejected sample is REPLACED by the median instead of being excluded. So the
            # test is exactly the one used by "sigma" above: deviation about the MEDIAN of
            # the current samples, in units of the unweighted sample standard deviation
            # (N-1) about their MEAN, both recomputed every pass on the already-modified
            # stack. (Established from Siril's documentation, not from rejection_float.c.)
            # Unlike "sigma" there is no zero-median guard: it is not known whether Siril
            # applies one to this mode, and the previous behaviour had none. Replacement
            # never lowers the sample count, so the 4-sample limit does not apply.
            if not np.any(can_attempt_rejection):
                break

            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'All-NaN slice')
                med_replace = np.nanmedian(np.where(active, stack, np.nan), axis=0)

            dev = stack - med_replace
            sigma = _active_sample_sd(stack, active) + 1e-12
            newly_bad = active & ((dev < -sigma_low * sigma) | (dev > sigma_high * sigma))
            newly_bad = newly_bad & can_attempt_rejection     # eligible pixels only
            if not newly_bad.any():
                break                                          # converged: nothing to replace
            stack = np.where(newly_bad, med_replace, stack)

        elif rejection == "percentile":
            # Percentile clipping matching Siril's rejection_float.c:177-190
            # Non-iterative, percentage-based thresholds relative to median
            
            # Compute median for zero-check (Siril lines 168-170)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'All-NaN slice')
                med_check = np.nanmedian(np.where(active, stack, np.nan), axis=0)
            
            # Siril lines 168-170: "if (median == 0.0) return 0"
            median_is_valid = (med_check != 0.0) & np.isfinite(med_check)
            can_reject_percentile = can_attempt_rejection & median_is_valid
            
            if not np.any(can_reject_percentile):
                break
            
            # Apply percentile clipping
            for frame_idx in range(n):
                if not active[frame_idx].any():
                    continue
                
                pixel_vals = stack[frame_idx]
                frame_active = active[frame_idx]
                
                # Vectorized percentile test
                low_reject = frame_active & ((med_check - pixel_vals) > med_check * sigma_low)
                high_reject = frame_active & ((pixel_vals - med_check) > med_check * sigma_high)
                
                # Only reject on eligible pixels
                low_reject = low_reject & can_reject_percentile
                high_reject = high_reject & can_reject_percentile
                
                active[frame_idx] = frame_active & ~low_reject & ~high_reject
            
            break  # Non-iterative
                
        elif rejection == "linearfit":
            # Iterates internally until a pass rejects nothing.
            active = _linearfit_reject_numpy(stack, active, sigma_low, sigma_high)
            break

        elif rejection == "gesdt":
            # Generalized ESD test matching Siril's rejection_float.c:318-355
            # Statistically rigorous outlier detection using Rosner's method
            
            if not SCIPY_AVAILABLE:
                raise RuntimeError("GESDT rejection requires scipy. Install with: pip install scipy")
            
            # Vectorized over all pixels of the strip (sort along the frame axis, then a
            # short loop over outlier candidates).
            active = _gesdt_reject_numpy(stack, active, gesdt_alpha, gesdt_max_frac)
            break  # Non-iterative

    # Final weighted statistics
    wa = np.where(active, w, 0.0)
    total_w = wa.sum(axis=0)
    total_w_safe = np.where(total_w > 0, total_w, 1.0)
    mean = (wa * stack).sum(axis=0) / total_w_safe
    m2 = (wa * (stack - mean) ** 2).sum(axis=0)
    count = active.sum(axis=0).astype(np.uint16)
    weight = total_w.astype(np.float32)      # exactly 0 wherever count == 0
    w2 = (wa * wa).sum(axis=0).astype(np.float32)    # sum of squared weights of the accepted samples

    return mean.astype(np.float32), m2.astype(np.float32), weight, w2, count


# --------------------------------------------------------------------------
# Optional PyTorch backend for compute_batch_stats
#
# All per-pixel work (sorting, medians, weighted sums, Winsorizing, linear-fit
# and GESDT rejection) is expressed as tensor ops over the frame axis, so every
# pixel of a strip is processed in parallel on a GPU (CUDA / Apple MPS) or on
# multi-core CPU tensors. PyTorch is located the same way amsp_LSF.py does it:
# sirilpy.TorchHelper().ensure_torch(), then `import torch`. If PyTorch is
# missing, disabled, or raises at runtime, the NumPy implementation is used.
# --------------------------------------------------------------------------

_TORCH_STATE = {"tried": False, "ok": False, "device": None, "pref": "auto",
                "max_samples": None, "oom_splits": 0, "fallbacks": 0, "seen": set(),
                "precision": "auto"}
_ESD_TABLE_CACHE: dict = {}
_LOGGER = [None]


def set_logger(fn) -> None:
    """Route this module's warnings to `fn(text)` (e.g. siril.log); print() if None."""
    _LOGGER[0] = fn


def _say(msg: str) -> None:
    fn = _LOGGER[0]
    try:
        (fn or print)(msg)
    except Exception:
        print(msg)


def bytes_per_sample(rejection: str) -> int:
    """Working-memory estimate per (frame, pixel) sample for the strip-size budget in
    lsf_stack.py. torch keeps sort indices and masks alive (96). NumPy peaks measured with
    tracemalloc (30 frames, float32 input): sigma/mad/median_sigma/percentile/gesdt 37-47
    (gesdt works in column blocks, so its temporaries do not grow with the strip),
    winsorized 69, linearfit 92; the figures below add about 15% headroom."""
    if _resolve_torch():
        return 96
    return {"linearfit": 106, "winsorized": 80}.get(rejection, 48)


def torch_summary() -> "str | None":
    """One-line recap of OOM splits / NumPy fallbacks, or None if there were none."""
    st = _TORCH_STATE
    if not (st["oom_splits"] or st["fallbacks"]):
        return None
    return (f"PyTorch backend: {st['oom_splits']} out-of-memory strip split(s)"
            + (f" (strips now capped at {st['max_samples']:,} samples)" if st["max_samples"] else "")
            + f", {st['fallbacks']} strip(s) fell back to NumPy.")


def configure_torch(device: str = "auto", precision: str = "auto") -> str:
    """Select the compute backend. device: auto | cpu | cuda | mps | off.
    precision: auto | float32 | float64 -- the dtype of the rejection work (sort, median,
    thresholds). auto = float32 on CUDA (consumer GPUs run float64 at 1/32 to 1/64 of
    the float32 rate) and MPS (no float64 there), float64 on CPU tensors. The final
    weighted sums are float64 whenever the device has it, whatever this is set to.
    Returns a short human-readable description of what will be used."""
    precision = (precision or "auto").lower()
    if precision not in ("auto", "float32", "float64"):
        raise ValueError(f"unknown torch precision '{precision}' (auto | float32 | float64)")
    _TORCH_STATE.update(tried=False, ok=False, device=None,
                        pref=(device or "auto").lower(),
                        max_samples=None, oom_splits=0, fallbacks=0, seen=set(),
                        precision=precision)
    _resolve_torch()
    return torch_description()


def _resolve_torch() -> bool:
    st = _TORCH_STATE
    if st["tried"]:
        return st["ok"]
    st["tried"] = True
    pref = st["pref"]
    if pref == "off":
        return False
    try:
        try:
            import sirilpy as _s
            _s.TorchHelper().ensure_torch()
        except Exception:
            pass   # not running inside Siril / helper unavailable: plain import below
        import torch
    except Exception:
        return False
    dev = None
    try:
        if pref in ("auto", "cuda") and torch.cuda.is_available():
            dev = torch.device("cuda")
        elif (pref in ("auto", "mps") and hasattr(torch.backends, "mps")
              and torch.backends.mps.is_available()):
            dev = torch.device("mps")
        elif pref in ("auto", "cpu"):
            dev = torch.device("cpu")
    except Exception:
        dev = None
    if dev is None:
        return False
    st["device"], st["ok"] = dev, True
    return True


def _torch_dtypes(torch, dev):
    """(compute dtype, accumulation dtype) for this device and the precision setting.
    MPS has no float64, so both are float32 there."""
    if dev.type == "mps":
        return torch.float32, torch.float32
    pref = _TORCH_STATE.get("precision", "auto")
    if pref == "float32" or (pref == "auto" and dev.type == "cuda"):
        return torch.float32, torch.float64
    return torch.float64, torch.float64


def torch_description() -> str:
    if not _resolve_torch():
        return "NumPy (CPU) -- PyTorch off or unavailable"
    import torch
    dev = _TORCH_STATE["device"]
    dt, acc = _torch_dtypes(torch, dev)
    prec = (f"{str(dt).replace('torch.', '')} compute, {str(acc).replace('torch.', '')} sums"
            if dt != acc else f"{str(dt).replace('torch.', '')}")
    if dev.type == "cuda":
        extra = f" ({torch.cuda.get_device_name(0)}, {prec})"
    elif dev.type == "mps":
        extra = f" ({prec}; MPS has no float64)"
    else:
        extra = f" ({prec})"
    return f"PyTorch {torch.__version__} on {dev.type}{extra}"


def _t_masked_median(torch, x, active):
    """Median over dim 0 of the active samples (even counts average the two
    middle values, like np.nanmedian). NaN where nothing is active."""
    n = x.shape[0]
    xs = torch.where(active, x, torch.full_like(x, float("inf"))).sort(dim=0).values
    k = active.sum(0)
    lo = torch.div((k - 1).clamp_min(0), 2, rounding_mode="floor").unsqueeze(0)
    hi = torch.div(k, 2, rounding_mode="floor").clamp_max(n - 1).unsqueeze(0)
    med = (xs.gather(0, lo).squeeze(0) + xs.gather(0, hi).squeeze(0)) * 0.5
    return torch.where(k > 0, med, torch.full_like(med, float("nan")))


def _t_final_stats(torch, x, active, w64, acc, max_elems=4_000_000):
    """Final weighted mean / M2 / total weight / sum of squared weights / count in the accumulation dtype `acc`
    (float64 where the device has it), done in frame chunks so the float64 copy of the
    strip never exists at once: only chunk x P samples are widened at a time. x is zero
    where inactive; w64 is (n, 1). Two passes (mean, then M2 about it), as in the numpy
    path."""
    n, P = x.shape
    chunk = max(1, min(n, max_elems // max(P, 1)))
    zero = torch.zeros((), device=x.device, dtype=acc)
    tw = torch.zeros(P, device=x.device, dtype=acc)
    tw2 = torch.zeros(P, device=x.device, dtype=acc)
    sx = torch.zeros(P, device=x.device, dtype=acc)
    for i in range(0, n, chunk):
        wc = torch.where(active[i:i + chunk], w64[i:i + chunk], zero)
        tw += wc.sum(0)
        tw2 += (wc * wc).sum(0)
        sx += (wc * x[i:i + chunk].to(acc)).sum(0)
    mean = sx / torch.where(tw > 0, tw, torch.ones((), device=x.device, dtype=acc))
    m2 = torch.zeros(P, device=x.device, dtype=acc)
    for i in range(0, n, chunk):
        wc = torch.where(active[i:i + chunk], w64[i:i + chunk], zero)
        dv = x[i:i + chunk].to(acc) - mean
        m2 += (wc * dv * dv).sum(0)
    return mean, m2, tw, tw2, active.sum(0)


def _t_active_sample_sd(torch, x, active, dt):
    """Torch twin of _active_sample_sd: unweighted sample SD (N-1) about the mean of
    the active samples, per pixel (column of the flat (N, P) tensors). `x` is already
    zero where inactive."""
    act_f = active.to(dt)
    kf = act_f.sum(0)
    m_u = (x * act_f).sum(0) / kf.clamp_min(1)
    dm = (x - m_u) * act_f
    return ((dm * dm).sum(0) / (kf - 1).clamp_min(1)).sqrt()


def _t_limit_bad(torch, bad, k, dev, sigma, sigma_low, sigma_high, min_keep=MIN_KEEP):
    """Torch twin of _limit_bad_numpy: per pixel (column of the flat (N, P)
    tensors) allow at most max(k - min_keep, 0) rejections, keeping the most
    extreme candidates. Only columns over budget are ranked/scored."""
    n = bad.shape[0]
    budget = (k - min_keep).clamp_min(0)
    over = bad.sum(0) > budget
    if not bool(over.any()):
        return bad
    cols = over.nonzero().squeeze(1)
    sub_bad = bad[:, cols]                                              # (N, M)
    d = dev[:, cols]
    thr = torch.where(d < 0, torch.full_like(d, sigma_low), torch.full_like(d, sigma_high)) * sigma[cols].unsqueeze(0)
    score = d.abs() / thr
    sub_s = torch.where(sub_bad, score, torch.full_like(score, float("-inf")))
    order = sub_s.argsort(dim=0, descending=True)                       # most extreme first
    rank = torch.empty_like(order)
    rank.scatter_(0, order, torch.arange(n, device=bad.device).view(n, 1).expand(n, order.shape[1]).contiguous())
    out = bad.clone()
    out[:, cols] = sub_bad & (rank < budget[cols].unsqueeze(0))
    return out


def _t_finish_reject(torch, active, newly_bad, dev, sigma, sigma_low, sigma_high):
    """Apply the 'never leave fewer than MIN_KEEP samples' rule (Siril: stop
    once N - r <= 4) per pixel, rejecting the most extreme candidates first
    when the limit binds. Returns the new active mask."""
    return active & ~_t_limit_bad(torch, newly_bad, active.sum(0), dev, sigma, sigma_low, sigma_high)


def _esd_crit_table(n: int, alpha: float) -> "np.ndarray":
    """Critical values indexed by sample size (0..n); inf where undefined."""
    key = (n, float(alpha))
    tab = _ESD_TABLE_CACHE.get(key)
    if tab is None:
        tab = np.full(n + 1, np.inf)
        s = np.arange(3, n + 1, dtype=np.float64)
        if s.size:
            with np.errstate(all="ignore"):
                tv = t_dist.ppf(1 - alpha / (2 * s), s - 2)
                vals = (s - 1) * tv / (np.sqrt(s) * np.sqrt(s - 2 + tv ** 2))
            tab[3:] = np.nan_to_num(vals, nan=np.inf)
        _ESD_TABLE_CACHE[key] = tab
    return tab


def compute_batch_stats_torch(stack, frame_weights, rejection="sigma",
                              sigma_low=3.0, sigma_high=3.0, max_iters=3,
                              gesdt_alpha=0.05, gesdt_max_frac=0.3):
    """PyTorch twin of compute_batch_stats_numpy (same inputs/outputs)."""
    import torch
    dev = _TORCH_STATE["device"]
    dt, acc = _torch_dtypes(torch, dev)     # rejection work in dt, final sums in acc
    eps_dt = torch.finfo(dt).eps
    tiny_dt = torch.finfo(dt).tiny          # 1e-300 underflows to 0 in float32
    n = stack.shape[0]
    out_shape = stack.shape[1:]

    flat = np.ascontiguousarray(stack).reshape(n, -1)
    P = flat.shape[1]
    x = torch.from_numpy(flat).to(device=dev, dtype=dt)
    w64 = torch.from_numpy(np.asarray(frame_weights, dtype=np.float64)).to(device=dev, dtype=acc).view(n, 1)
    w = w64.to(dt)
    zero = torch.zeros((), device=dev, dtype=dt)
    one = torch.ones((), device=dev, dtype=dt)
    idx = torch.arange(n, device=dev).view(n, 1)

    active = torch.isfinite(x)          # NaN = no data at this pixel
    x = torch.where(active, x, zero)
    can_attempt = active.sum(0) > 1

    # `iters` counts rejection passes; final statistics are recomputed after the loop.
    iters = 0 if (rejection == "none" or n < 3) else max_iters

    for it in range(iters):
        if rejection in ("mad", "sigma"):
            med = _t_masked_median(torch, x, active)
            can = can_attempt & (med != 0) & torch.isfinite(med)
            if not bool(can.any()):
                break
            dev_ = x - med
            if rejection == "mad":
                sigma = _t_masked_median(torch, (x - med).abs(), active) * 1.4826 + 1e-12
            else:
                # Siril's estimator: unweighted sample SD (N-1) about the mean of the
                # active samples; the clip test stays about the median (see numpy path).
                sigma = _t_active_sample_sd(torch, x, active, dt) + 1e-12
            bad = active & ((dev_ < -sigma_low * sigma) | (dev_ > sigma_high * sigma)) & can
            still = _t_finish_reject(torch, active, bad, dev_, sigma, sigma_low, sigma_high)
            if torch.equal(still, active):
                break
            active = still

        elif rejection == "winsorized":
            if not bool(can_attempt.any()):
                break
            medianf = _t_masked_median(torch, x, active)
            can_w = can_attempt & torch.isfinite(medianf)
            if not bool(can_w.any()):
                break
            medianf = torch.nan_to_num(medianf, nan=0.0)
            wa = torch.where(active, w, zero)            # weighted mean of the active samples
            tw_ = wa.sum(0)                              # (only this method needs it)
            tws = torch.where(tw_ > 0, tw_, one)
            mean = (wa * x).sum(0) / tws
            dv = x - mean
            sigma_ws = ((wa * dv * dv).sum(0) / tws).clamp_min(0).sqrt() + 1e-12
            w_stack = x.clone()
            for _ in range(20):
                prev = sigma_ws
                m0 = (medianf - 1.5 * sigma_ws).unsqueeze(0)
                m1 = (medianf + 1.5 * sigma_ws).unsqueeze(0)
                w_stack = torch.minimum(torch.maximum(w_stack, m0), m1)
                dvw = w_stack - mean
                sigma_ws = 1.134 * ((wa * dvw * dvw).sum(0) / tws).clamp_min(0).sqrt() + 1e-12
                rel = (sigma_ws - prev).abs() / (prev + 1e-12)
                if bool(((rel < 0.0005) | ~can_w).all()):
                    break
            dev_ = x - medianf
            bad = (active & ((dev_ < -sigma_low * sigma_ws) | (dev_ > sigma_high * sigma_ws))) & can_w
            still = _t_finish_reject(torch, active, bad, dev_, sigma_ws, sigma_low, sigma_high)
            if torch.equal(still, active):
                break
            active = still

        elif rejection == "median_sigma":
            # Same test as "sigma" (median-centred, unweighted N-1 SD about the mean), with
            # rejected samples replaced by the median instead of excluded. See the numpy path.
            if not bool(can_attempt.any()):
                break
            med = _t_masked_median(torch, x, active)
            dev_ = x - med
            sigma = _t_active_sample_sd(torch, x, active, dt) + 1e-12
            bad = active & ((dev_ < -sigma_low * sigma) | (dev_ > sigma_high * sigma)) & can_attempt
            if not bool(bad.any()):
                break                                    # converged: nothing to replace
            x = torch.where(bad, med.unsqueeze(0), x)    # replace, don't exclude

        elif rejection == "percentile":
            med = _t_masked_median(torch, x, active)
            can = can_attempt & (med != 0) & torch.isfinite(med)
            if not bool(can.any()):
                break
            m = med.unsqueeze(0)
            low = active & ((m - x) > m * sigma_low)
            high = active & ((x - m) > m * sigma_high)
            active = active & ~((low | high) & can.unsqueeze(0))
            break   # non-iterative

        elif rejection == "linearfit":
            # Batched version of Siril's per-pixel sorted linear-fit clipping.
            for _ in range(max(n, 1)):
                xs, order = torch.where(active, x, torch.full_like(x, float("inf"))).sort(0)
                k = active.sum(0)
                valid = idx < k.unsqueeze(0)
                ys = torch.where(valid, xs, zero)
                kf = k.to(dt).clamp_min(1)
                m_x = (kf - 1) / 2
                m_y = ys.sum(0) / kf
                di = idx.to(dt) - m_x.unsqueeze(0)
                dy = torch.where(valid, ys - m_y.unsqueeze(0), zero)
                var = torch.where(valid, di * di, zero).sum(0) / kf
                slope = ((di * dy).sum(0) / kf) / var.clamp_min(tiny_dt)
                icpt = m_y - m_x * slope
                expected = icpt.unsqueeze(0) + slope.unsqueeze(0) * idx.to(dt)
                sig = torch.where(valid, (ys - expected).abs(), zero).sum(0) / kf
                low = (expected - ys) > sig.unsqueeze(0) * sigma_low
                high = (ys - expected) > sig.unsqueeze(0) * sigma_high
                cand = valid & (low | high) & (k > 4).unsqueeze(0)
                # keep >= MIN_KEEP samples; most extreme candidates go first when over budget
                allowed = _t_limit_bad(torch, cand, k, ys - expected, sig, sigma_low, sigma_high)
                if not bool(allowed.any()):
                    break
                rej = torch.zeros((n, P), dtype=torch.long, device=dev)
                rej.scatter_(0, order, allowed.to(torch.long))
                active = active & (rej == 0)
            break

        elif rejection == "gesdt":
            if not SCIPY_AVAILABLE:
                raise RuntimeError("GESDT rejection requires scipy. Install with: pip install scipy")
            k = active.sum(0)
            # Outlier budget per pixel, computed on the host in float64 (floor(k * frac) must
            # not depend on the compute dtype, and MPS has no float64); only P integers move.
            k_np = k.cpu().numpy().astype(np.float64)
            mo = np.maximum(np.floor(k_np * gesdt_max_frac), 1.0)
            mo = np.where(mo >= k_np, k_np - 1.0, mo)
            mo = np.where(k_np > 1, mo, 0.0).astype(np.int64)
            max_out = torch.from_numpy(mo).to(dev)
            T = int(mo.max()) if mo.size else 0
            if T > 0:
                table = torch.from_numpy(_esd_crit_table(n, gesdt_alpha)).to(device=dev, dtype=dt)
                xs, order = torch.where(active, x, torch.full_like(x, float("inf"))).sort(0)
                valid = idx < k.unsqueeze(0)
                ys = torch.where(valid, xs, zero)
                lo = torch.zeros(P, dtype=torch.long, device=dev)
                hi = (k - 1).clamp_min(0)
                outs, poss = [], []
                for t in range(T):
                    vt = max_out > t
                    size = hi - lo + 1
                    maskw = (idx >= lo.unsqueeze(0)) & (idx <= hi.unsqueeze(0))
                    sz = size.clamp_min(1).to(dt)
                    mean_w = torch.where(maskw, ys, zero).sum(0) / sz
                    dv = torch.where(maskw, ys - mean_w.unsqueeze(0), zero)
                    std = ((dv * dv).sum(0) / (sz - 1).clamp_min(1)).sqrt()
                    y_lo = ys.gather(0, lo.unsqueeze(0)).squeeze(0)
                    y_hi = ys.gather(0, hi.unsqueeze(0)).squeeze(0)
                    dlow, dhigh = mean_w - y_lo, y_hi - mean_w
                    use_high = dhigh > dlow
                    # SD below a few ulps of the mean is rounding noise, not scatter (a constant
                    # pixel would otherwise get a huge random G); same idea as the numpy path.
                    good = std > (16.0 * eps_dt) * mean_w.abs() + tiny_dt
                    G = torch.where(good, torch.where(use_high, dhigh, dlow) / torch.where(good, std, one), zero)
                    outs.append((G > table[size.clamp(0, n)]) & vt)
                    poss.append(torch.where(use_high, hi, lo))
                    lo = lo + (vt & ~use_high).to(torch.long)
                    hi = hi - (vt & use_high).to(torch.long)
                is_out = torch.stack(outs)                       # (T, P)
                pos = torch.stack(poss)
                # Rosner: everything up to the LAST significant candidate is an outlier.
                marked = is_out.flip(0).to(torch.long).cumsum(0).flip(0) > 0
                marked = marked & (torch.arange(T, device=dev).view(T, 1) < max_out.unsqueeze(0))
                rej_sorted = torch.zeros((n, P), dtype=torch.long, device=dev)
                rej_sorted.scatter_add_(0, pos, marked.to(torch.long))
                rej = torch.zeros((n, P), dtype=torch.long, device=dev)
                rej.scatter_(0, order, rej_sorted)
                active = active & (rej == 0)
            break

    # Final weighted statistics, in the accumulation dtype (float64 on CUDA/CPU).
    mean, m2, tw, tw2, count = _t_final_stats(torch, x, active, w64, acc)

    def _np(t, dtype):
        return t.cpu().numpy().astype(dtype, copy=False).reshape(out_shape)

    return (_np(mean, np.float32), _np(m2, np.float32), _np(tw, np.float32), _np(tw2, np.float32),
            _np(count.clamp_max(65535), np.uint16))


def _is_oom(exc: BaseException) -> bool:
    """True for an out-of-memory failure from CUDA, MPS or the CPU allocator."""
    if isinstance(exc, MemoryError):
        return True
    text = (type(exc).__name__ + " " + str(exc)).lower()
    return any(s in text for s in ("out of memory", "outofmemory", "can't allocate memory",
                                   "cannot allocate memory", "insufficient memory"))


def _torch_free_cache() -> None:
    try:
        import torch
        dev = _TORCH_STATE["device"]
        if dev is not None and dev.type == "cuda":
            torch.cuda.empty_cache()
        elif dev is not None and dev.type == "mps":
            torch.mps.empty_cache()
    except Exception:
        pass


def _torch_stats_retry(stack, frame_weights, rejection, sigma_low, sigma_high, max_iters,
                       gesdt_alpha, gesdt_max_frac):
    """compute_batch_stats_torch with out-of-memory recovery: on OOM the strip is split in
    half along its longest spatial axis (pixels are independent, so the halves are exact)
    and each half is retried, recursively. The size that failed is remembered
    (`max_samples`), so later strips are split up front instead of failing again.
    Anything that is not an OOM, or an OOM on a single pixel, is re-raised."""
    st = _TORCH_STATE
    args = (frame_weights, rejection, sigma_low, sigma_high, max_iters, gesdt_alpha, gesdt_max_frac)
    spatial = [ax for ax in range(1, stack.ndim) if stack.shape[ax] > 1]

    def split():
        ax = max(spatial, key=lambda a: stack.shape[a])
        mid = stack.shape[ax] // 2
        pre = (slice(None),) * ax
        a = _torch_stats_retry(stack[pre + (slice(0, mid),)], *args)
        b = _torch_stats_retry(stack[pre + (slice(mid, None),)], *args)
        return tuple(np.concatenate([u, v], axis=ax - 1) for u, v in zip(a, b))

    cap = st["max_samples"]
    if spatial and cap and stack.size > cap:
        return split()
    try:
        return compute_batch_stats_torch(stack, *args)
    except Exception as exc:
        if not (spatial and _is_oom(exc)):
            raise
        _torch_free_cache()
        new_cap = max(1, (min(cap, stack.size) if cap else stack.size) // 2)
        st["max_samples"] = new_cap
        st["oom_splits"] += 1
        if ("oom", new_cap) not in st["seen"]:
            st["seen"].add(("oom", new_cap))
            _say(f"[lsf] PyTorch out of memory on {stack.size:,} samples; splitting strips "
                 f"(cap now {new_cap:,} samples). Lower --mem-budget-gb to avoid the retries.")
        return split()


def compute_batch_stats(stack, frame_weights, rejection="sigma",
                        sigma_low=3.0, sigma_high=3.0, max_iters=3,
                        gesdt_alpha=0.05, gesdt_max_frac=0.3):
    """Dispatcher: PyTorch backend when available/enabled, NumPy otherwise.
    A torch out-of-memory error is handled by splitting the strip and retrying on the GPU
    (see _torch_stats_retry). Any other torch failure falls back to NumPy for that call,
    and says so every time the message is new: for gesdt and linearfit the NumPy path is
    far slower, so a silent fallback can turn a minutes-long job into a long one."""
    if rejection not in SUPPORTED_REJECTION:
        raise ValueError(f"unsupported rejection method '{rejection}' (supported: {SUPPORTED_REJECTION})")
    if _resolve_torch():
        try:
            return _torch_stats_retry(stack, frame_weights, rejection,
                                      sigma_low, sigma_high, max_iters,
                                      gesdt_alpha, gesdt_max_frac)
        except Exception as exc:
            st = _TORCH_STATE
            st["fallbacks"] += 1
            key = ("fail", type(exc).__name__, str(exc)[:200])
            if key not in st["seen"]:
                st["seen"].add(key)
                _say(f"[lsf] WARNING: PyTorch backend failed ({type(exc).__name__}: {exc}); this strip "
                     f"falls back to NumPy, which is much slower for {rejection}.")
            _torch_free_cache()
    return compute_batch_stats_numpy(stack, frame_weights, rejection,
                                     sigma_low, sigma_high, max_iters,
                                     gesdt_alpha, gesdt_max_frac)


def combine_weighted_stats(mean_a, m2_a, w_a, w2_a, c_a, mean_b, m2_b, w_b, w2_b, c_b):
    """Parallel/pairwise merge of two weighted Welford accumulators
    (Chan et al.'s pairwise update, generalized to weights). This is the
    "numerically stable incremental update" referenced in Section 9 of the
    LSF spec, applied at the batch level rather than frame-by-frame so a
    whole new batch can be folded into history in one pass per pixel.

    The sum of squared weights (LSF_W2) is additive: every accepted sample
    contributes its own squared weight, so the merged W2 is simply W2_a + W2_b.
    Returns (mean, m2, weight, w2, count).
    """
    w_a = w_a.astype(np.float64)
    w_b = w_b.astype(np.float64)
    combined_w = w_a + w_b
    combined_w2 = w2_a.astype(np.float64) + w2_b.astype(np.float64)
    safe_w = np.where(combined_w > 0, combined_w, 1.0)
    delta = mean_b.astype(np.float64) - mean_a.astype(np.float64)
    combined_mean = mean_a + delta * (w_b / safe_w)
    combined_m2 = m2_a + m2_b + (delta ** 2) * (w_a * w_b / safe_w)
    combined_count = (c_a.astype(np.uint32) + c_b.astype(np.uint32))
    combined_count = np.minimum(combined_count, 65535).astype(np.uint16)
    return (
        combined_mean.astype(np.float32),
        combined_m2.astype(np.float32),
        combined_w.astype(np.float32),
        combined_w2.astype(np.float32),
        combined_count,
    )


def compute_frame_weight(method: str, noise: Optional[float], nbstars: Optional[int],
                          wfwhm: Optional[float], scale: float = 1.0) -> float:
    """Per-frame, per-channel stacking weight.

    `scale` is the ingestion-time normalization scale applied to this channel
    (Section 8.3). Siril's `bgnoise` is measured on the frame BEFORE
    normalization, but the accumulator holds the frame AFTER `pixel*scale +
    offset`, so the noise the weight must describe is `noise * |scale|`
    (the offset does not change it). Only the "noise" method is in flux
    units; nbstars and wfwhm are scale-free and ignore `scale`.
    """
    if method == "none":
        return 1.0
    if method == "noise":
        if not noise or noise <= 0:
            return 1.0
        eff = noise
        if scale is not None and np.isfinite(scale) and scale != 0:
            eff = noise * abs(float(scale))
        return 1.0 / (eff * eff)
    if method == "nbstars":
        if not nbstars or nbstars <= 0:
            return 1.0
        return float(nbstars)
    if method == "wfwhm":
        if not wfwhm or wfwhm <= 0:
            return 1.0
        return 1.0 / (wfwhm * wfwhm)
    raise ValueError(f"unknown weighting method '{method}'")


# --------------------------------------------------------------------------
# Derived statistics and per-pixel uncertainty maps (FITS profile, Section 10).
#
# With W = LSF_WEIGHT, M2 = LSF_M2 and W2 = LSF_W2 (sum of squared weights) at
# a pixel and plane:
#
#   n_eff       = W^2 / W2             effective sample size
#   s^2         = M2 / (W - W2/W)      unbiased variance of ONE sample about the mean
#                                      (undefined with fewer than two contributing samples)
#   var(mean)   = s^2 * W2 / W^2       variance of the stacked mean: the noise left in
#                                      the exported master
#   var(new)    = s^2 * (1 + W2/W^2)   expected variance of a NEW sample about the
#                                      accumulated mean (for a history rejection test)
#
# These assume the samples share a common underlying variance. Check against the
# textbook case of N equal weights (W = W2 = N): n_eff = N, s^2 = M2/(N-1),
# var(mean) = s^2/N, var(new) = s^2 (1 + 1/N).
#
# s^2 and var(mean) differ by about a factor n_eff (so sqrt(N) in noise): the
# sample maps describe the scatter of one sub-frame, the mean maps the
# uncertainty of the exported master.
#
# W2 is stored exactly, so n_eff is exact. It replaces the v0.2-era
# approximation count * weighting-efficiency, which was only right where every
# frame covered the pixel and nothing was rejected, and which the uint16
# saturation of the count plane could corrupt.
# --------------------------------------------------------------------------

ERROR_MAP_KINDS = ("sample_variance", "sample_noise", "mean_variance", "mean_noise", "n_eff")

# A pixel needs n_eff clearly above 1 for a variance to be measurable. A single sample has
# W - W2/W == 0 exactly, which float32 storage turns into rounding noise; this margin keeps
# that noise out (two samples whose weights differ by 1:20000 are "one sample" in practice).
NEFF_MIN = 1.0 + 1e-4


def _f64(a) -> np.ndarray:
    return np.asarray(a, dtype=np.float64)


def effective_sample_size(weight, w2) -> np.ndarray:
    """n_eff = W^2 / W2 (float64). 0 where the pixel has no data (W or W2 is 0)."""
    W, W2 = _f64(weight), _f64(w2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where((W > 0) & (W2 > 0), W * W / W2, 0.0)


def unbiased_variance(m2, weight, w2) -> np.ndarray:
    """s^2 = M2 / (W - W2/W) (float64): the unbiased variance of a single sample about the
    accumulated mean. NaN where it is undefined (n_eff <= 1, i.e. fewer than two samples,
    or no data)."""
    W, W2, M = _f64(weight), _f64(w2), _f64(m2)
    n_eff = effective_sample_size(W, W2)
    with np.errstate(divide="ignore", invalid="ignore"):
        s2 = M / (W - W2 / W)
    return np.where(n_eff > NEFF_MIN, s2, np.nan)


def variance_of_mean(m2, weight, w2) -> np.ndarray:
    """s^2 * W2 / W^2 (float64): variance of the stacked mean. NaN where s^2 is undefined."""
    W, W2 = _f64(weight), _f64(w2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return unbiased_variance(m2, W, W2) * W2 / (W * W)


def expected_new_sample_variance(m2, weight, w2) -> np.ndarray:
    """s^2 * (1 + W2/W^2) (float64): the variance to expect of a NEW sample about the
    accumulated mean, i.e. the scale for a z-score test of incoming frames against the
    stored history (the history rejection test of Section 10; its order relative to the
    batch-local rejection and its thresholds are still open, Section 13). NaN where s^2
    is undefined."""
    W, W2 = _f64(weight), _f64(w2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return unbiased_variance(m2, W, W2) * (1.0 + W2 / (W * W))


def error_map(kind: str, m2, weight, w2, rows: int = 256) -> np.ndarray:
    """Per-pixel uncertainty / confidence map as float32, computed in row bands (workers
    write disjoint slices of the output; no band results are kept alive).

    kind: "sample_variance" | "sample_noise" -> s^2 (and its sqrt): scatter of ONE frame
          "mean_variance"   | "mean_noise"   -> var(mean) (and its sqrt): noise of the stack
          "n_eff"                            -> effective sample size W^2 / W2
                                                (an effective-integration map, in frames)

    Conventions: pixels with no data (W == 0) are 0, matching the master. For the
    variance and noise maps a pixel with a single contributing sample has no measurable
    scatter, so it is NaN rather than a misleading 0. n_eff is 1 there.
    """
    if kind not in ERROR_MAP_KINDS:
        raise ValueError(f"unknown error map kind '{kind}' (supported: {ERROR_MAP_KINDS})")
    shape = m2.shape
    out = np.empty(shape, dtype=np.float32)
    take_sqrt = kind.endswith("noise")

    def band(r0):
        r1 = min(r0 + rows, shape[0])
        w = np.asarray(weight[r0:r1], dtype=np.float64)
        ww = np.asarray(w2[r0:r1], dtype=np.float64)
        if kind == "n_eff":
            out[r0:r1] = effective_sample_size(w, ww)
            return
        m = np.asarray(m2[r0:r1], dtype=np.float64)
        v = variance_of_mean(m, w, ww) if kind.startswith("mean") else unbiased_variance(m, w, ww)
        v = np.where(w > 0, v, 0.0)                       # no data -> 0 (NaN stays NaN: single sample)
        v = np.maximum(v, 0.0)                            # M2 rounding can dip a hair below 0
        out[r0:r1] = np.sqrt(v) if take_sqrt else v

    starts = list(range(0, shape[0], rows))
    n_workers = min(os.cpu_count() or 1, len(starts))
    if n_workers <= 1:
        for r0 in starts:
            band(r0)
    else:
        with ThreadPoolExecutor(max_workers=n_workers) as ex:
            list(ex.map(band, starts))
    return out


# --------------------------------------------------------------------------
# Small JSON convenience wrappers used by both scripts
# --------------------------------------------------------------------------

def dumps(obj) -> str:
    return json.dumps(obj, indent=2, default=str)


def loads(text: str):
    return json.loads(text)


# --------------------------------------------------------------------------
# Output normalization (Siril norm_to_0_1_range) stored as an invertible
# affine map, so the Stack plane can hold the Siril-matching image while
# SumSq/WeightSum stay on the raw scale. Parameters live in the
# "LSF:OutputNorm" property and the LSFNORM/LSFNMIN/LSFNMAX keywords.
# --------------------------------------------------------------------------

def apply_output_norm(mean: np.ndarray, weight: np.ndarray, rows: int = 256):
    """In-place min-max normalization of `mean` to [0,1] over covered pixels
    (weight > 0). min/max are taken over covered, non-zero pixels across all
    channels (as Siril does). Uncovered pixels stay 0. Returns (mini, maxi),
    or None if there is nothing to normalize."""
    mini, maxi = np.inf, -np.inf
    for r0 in range(0, mean.shape[0], rows):
        m = mean[r0:r0 + rows]
        v = (weight[r0:r0 + rows] > 0) & (m != 0) & np.isfinite(m)
        if v.any():
            vals = m[v]
            mini = min(mini, float(vals.min()))
            maxi = max(maxi, float(vals.max()))
    if not (np.isfinite(mini) and np.isfinite(maxi)) or maxi <= mini:
        return None
    rng = maxi - mini
    for r0 in range(0, mean.shape[0], rows):
        m = mean[r0:r0 + rows]           # view: writes go back into `mean`
        v = weight[r0:r0 + rows] > 0
        m[v] = ((m[v].astype(np.float64) - mini) / rng).astype(np.float32)
        m[~v] = 0.0
    return mini, maxi


def invert_output_norm(disp: np.ndarray, weight: np.ndarray, mini: float, maxi: float) -> np.ndarray:
    """Recover the raw weighted mean from a normalized Stack strip."""
    raw = disp.astype(np.float64) * (maxi - mini) + mini
    return np.where(weight > 0, raw, 0.0)


LSF_IMAGES = ("Stack", "WeightSum", "SumSq", "SumW2", "ContributionCount")
LSF_PROPERTIES = ("LSF:FrameLedger", "LSF:CanonicalGrid", "LSF:Settings", "LSF:OutputNorm")


def read_lsf(path: str) -> dict:
    """Read an LSF FITS file: {"images", "properties", "meta", "format"}. Planes are memory-mapped
    where the file allows it (big-endian; see native()). No validation: call validate_lsf_file()
    before merging into the file."""
    import lsf_fits
    if not lsf_fits.is_fits_file(path):
        raise ValueError(f"'{path}' is not a FITS file")
    return lsf_fits.read_lsf_fits(path)


def validate_lsf_file(path: str) -> None:
    """Section 8 pre-merge validation. Raises lsf_fits.LsfValidationError naming the failed check."""
    import lsf_fits
    lsf_fits.validate_lsf_fits(path)


def native(arr) -> np.ndarray:
    """C-contiguous native-byte-order copy (FITS planes are big-endian memory maps)."""
    a = np.asarray(arr)
    return np.ascontiguousarray(a, dtype=a.dtype.newbyteorder("="))


def write_lsf(path: str, images: list, properties: dict, compress=None, meta: "dict | None" = None,
              new_rows: "int | None" = None, row_order: "str | None" = None) -> dict:
    """Write an LSF FITS file. Returns {"bytes", "mode", "uid", "gen"}.

    compress   falsy | "gzip2" (lossless tile compression of the extension images)
    meta       read_lsf()["meta"] of the file being appended to: carries LSFUID, LSFGEN and everything
               that must survive a rewrite (C7). None when initializing.
    new_rows   number of trailing ledger rows added by this commit.
    """
    import lsf_fits
    return lsf_fits.write_lsf_fits(path, images, properties, compress=compress, meta=meta,
                                   new_rows=new_rows, row_order=row_order)


def check_lsf_file(data: dict, path: str, optional_images: tuple = ()) -> None:
    """Raise RuntimeError unless `data` (from read_lsf) carries every image and metadata property
    an LSF file has. `optional_images` names images the caller can do without (the exporter does
    not need SumW2 to export the master or the count map)."""
    missing = ([f"'{n}' image" for n in LSF_IMAGES
                if n not in optional_images and n not in data["images"]]
               + [p for p in LSF_PROPERTIES if p not in data["properties"]])
    if missing:
        hint = (" The LSF_W2 plane (sum of squared weights) is new in FITS profile v0.3; a file written "
                "before it existed can still be read as an image but cannot be merged into or used for "
                "noise maps." if "'SumW2' image" in missing else "")
        raise RuntimeError(f"'{path}' is not a valid LSF file (missing {', '.join(missing)}).{hint}")


def read_lsf_settings(path: str) -> dict:
    """Cheap peek at an existing LSF file's recorded processing options, for
    pre-filling a UI or inheriting defaults on --append. Returns
    {"settings": {...}, "n_frames": int}."""
    d = read_lsf(path)
    check_lsf_file(d, path)
    props = d["properties"]
    out = {
        "settings": loads(props["LSF:Settings"]),
        "n_frames": len(loads(props["LSF:FrameLedger"])),
    }
    del d
    return out