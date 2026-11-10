/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/* LSF frame weights (FITS profile v0.4, Section 10.3).
 * Siril's per-frame inputs and formulas (compute_noise_weights,
 * compute_wfwhm_weights and compute_nbstars_weights in median_and_mean.c),
 * but with the constants frozen from the first batch, so that weights stay
 * comparable across appends, and with a floor (weight.floor) where
 * Siril gives the worst frame 0. Siril's own functions are not changed. */

#include <float.h>
#include <limits.h>
#include <math.h>

#include "core/siril_log.h"
#include "lsf/lsf.h"
#include "lsf/lsf_weights.h"

double lsf_noise_raw_weight(double pscale, double bgnoise) {
	double d = pscale * pscale * bgnoise * bgnoise;
	return d > 0.0 ? 1.0 / d : 0.0;
}

static double wfwhm_formula(const lsf_settings *s, float wfwhm) {
	if (wfwhm <= 0.f)
		return 0.0;
	if (s->weight_model == LSF_WMODEL_INVERSE) {
		/* the limit of Siril's formula when fwhm_max is infinite: frames
		 * better than the first batch's best still weigh more */
		double r = s->fwhm_min / wfwhm;
		return MAX(r * r, s->weight_floor);
	}
	if (s->fwhm_max <= s->fwhm_min)
		return 1.0;
	double f = wfwhm;
	if (f < s->fwhm_min) f = s->fwhm_min;
	if (f > s->fwhm_max) f = s->fwhm_max;
	double invfwhmax2 = 1.0 / (s->fwhm_max * s->fwhm_max);
	double invdenom = 1.0 / (1.0 / (s->fwhm_min * s->fwhm_min) - invfwhmax2);
	return MAX((1.0 / (f * f) - invfwhmax2) * invdenom, s->weight_floor);
}

static double nbstars_formula(const lsf_settings *s, int nstars) {
	if (s->weight_model == LSF_WMODEL_INVERSE) {
		/* unclamped: a frame with more stars than the first batch's most
		 * still weighs more */
		double r = s->stars_max > 0.0 ? nstars / s->stars_max : 1.0;
		return MAX(r * r, s->weight_floor);
	}
	if (s->stars_max <= s->stars_min)
		return 1.0;
	double n = nstars;
	if (n < s->stars_min) n = s->stars_min;
	if (n > s->stars_max) n = s->stars_max;
	double invdenom = 1.0 / (s->stars_max - s->stars_min);
	return MAX((n - s->stars_min) * (n - s->stars_min) * invdenom * invdenom, s->weight_floor);
}

/* freezes the weight constants from the first batch.
 * noise_raw is indexed [chan * nb_frames + i] */
int lsf_freeze_weights(lsf_settings *s, int nb_frames, int nchan, const double *noise_raw,
		const float *wfwhm, const int *nstars) {
	if (nb_frames <= 0)
		return 1;
	switch (s->weighting) {
	case NOISE_WEIGHT:
		for (int c = 0; c < nchan; c++) {
			double norm = 0.0;
			for (int i = 0; i < nb_frames; i++)
				norm += noise_raw[c * nb_frames + i];
			norm /= nb_frames;
			if (norm <= 0.0) {
				siril_log_error(_("LSF: cannot compute noise weights (no background noise estimate)\n"));
				return 1;
			}
			s->noise_norm[c] = norm;
		}
		break;
	case WFWHM_WEIGHT: {
		double fmin = DBL_MAX, fmax = -DBL_MAX;
		for (int i = 0; i < nb_frames; i++) {
			if (wfwhm[i] > 0.f && wfwhm[i] < fmin) fmin = wfwhm[i];
			if (wfwhm[i] > fmax) fmax = wfwhm[i];
		}
		if (fmin == DBL_MAX || fmax <= 0.0) {
			siril_log_error(_("Sequence does not have registration info, cannot use weighing by %s, aborting\n"), "wFWHM");
			return 1;
		}
		/* weight.range widens the span of the first batch toward the worse
		 * frames: 1 is Siril's rule */
		s->fwhm_min = fmin;
		s->fwhm_max = fmin + s->weight_range * (fmax - fmin);
		s->fwhm_norm = 1.0;
		double norm = 0.0;
		for (int i = 0; i < nb_frames; i++)
			norm += wfwhm_formula(s, wfwhm[i]);
		norm /= nb_frames;
		if (norm <= 0.0)
			return 1;
		s->fwhm_norm = norm;
		break;
	}
	case NBSTARS_WEIGHT: {
		int nmin = INT_MAX, nmax = 0;
		for (int i = 0; i < nb_frames; i++) {
			if (nstars[i] < nmin) nmin = nstars[i];
			if (nstars[i] > nmax) nmax = nstars[i];
		}
		s->stars_max = nmax;
		s->stars_min = MAX(0.0, nmax - s->weight_range * (nmax - nmin));
		s->stars_norm = 1.0;
		double norm = 0.0;
		for (int i = 0; i < nb_frames; i++)
			norm += nbstars_formula(s, nstars[i]);
		norm /= nb_frames;
		if (norm <= 0.0)
			return 1;
		s->stars_norm = norm;
		break;
	}
	case NBSTACK_WEIGHT:
		siril_log_error(_("LSF: weighting by number of stacked images is not supported\n"));
		return 1;
	default:
		break;
	}
	return 0;
}

/* The frames of a batch outside the frozen range of a range model: below
 * the first batch's best (they all get the top weight, whatever their
 * differences) and beyond its worst (they all get the weight of the worst,
 * or the floor). Returns FALSE when the scheme or the model has no range. */
gboolean lsf_weight_range_check(const lsf_settings *s, int nb, const float *wfwhm, const int *nstars,
		int *better, int *worse, int *at_floor) {
	*better = *worse = *at_floor = 0;
	if (s->weight_model != LSF_WMODEL_RANGE || (s->weighting != WFWHM_WEIGHT && s->weighting != NBSTARS_WEIGHT))
		return FALSE;
	for (int i = 0; i < nb; i++) {
		double raw;
		if (s->weighting == WFWHM_WEIGHT) {
			if (wfwhm[i] <= 0.f)
				continue;
			*better += wfwhm[i] < s->fwhm_min;
			*worse += wfwhm[i] > s->fwhm_max;
			raw = wfwhm_formula(s, wfwhm[i]);
		} else {
			*better += nstars[i] > s->stars_max;
			*worse += nstars[i] < s->stars_min;
			raw = nbstars_formula(s, nstars[i]);
		}
		*at_floor += raw <= s->weight_floor;
	}
	return TRUE;
}

/* Switches a range model to inverse (Section 10.3), from this append on.
 * The frames already merged keep their weights, which the accumulator holds.
 * The new norm is chosen so that, under inverse, they would weigh in total
 * what they weigh now: hist_wfwhm, hist_nstars and hist_weight are the
 * ledger's WFWHM, NSTARS and WEIGHT of the first channel, n rows. The new
 * frames are then on the same scale as the history on average, without
 * rescaling the accumulator. Returns 0, or 1 when the ledger cannot give
 * the norm. */
int lsf_switch_to_inverse(lsf_settings *s, int n, const float *hist_wfwhm, const int *hist_nstars,
		const double *hist_weight) {
	if (s->weighting != WFWHM_WEIGHT && s->weighting != NBSTARS_WEIGHT)
		return 1;
	int model = s->weight_model;
	s->weight_model = LSF_WMODEL_INVERSE;
	double sum_inv = 0.0, sum_w = 0.0;
	for (int i = 0; i < n; i++) {
		if (!(hist_weight[i] > 0.0))
			continue;
		double raw = s->weighting == WFWHM_WEIGHT ? (hist_wfwhm[i] > 0.f ? wfwhm_formula(s, hist_wfwhm[i]) : 0.0)
			: nbstars_formula(s, hist_nstars[i]);
		if (!(raw > 0.0))
			continue;
		sum_inv += raw;
		sum_w += hist_weight[i];
	}
	if (!(sum_inv > 0.0) || !(sum_w > 0.0)) {
		s->weight_model = model;
		return 1;
	}
	if (s->weighting == WFWHM_WEIGHT)
		s->fwhm_norm = sum_inv / sum_w;
	else s->stars_norm = sum_inv / sum_w;
	return 0;
}

double lsf_frame_weight(const lsf_settings *s, int chan, double noise_raw, float wfwhm, int nstars) {
	switch (s->weighting) {
	case NOISE_WEIGHT:
		return s->noise_norm[chan] > 0.0 ? noise_raw / s->noise_norm[chan] : 0.0;
	case WFWHM_WEIGHT:
		return wfwhm_formula(s, wfwhm) / s->fwhm_norm;
	case NBSTARS_WEIGHT:
		return nbstars_formula(s, nstars) / s->stars_norm;
	default:
		return 1.0;
	}
}
