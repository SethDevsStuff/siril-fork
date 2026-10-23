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
 * comparable across appends, and with a floor (LSF_WEIGHT_FLOOR) where
 * Siril gives the worst frame 0. Siril's own functions are not changed. */

#include <math.h>
#include <float.h>
#include <limits.h>

#include "core/siril_log.h"
#include "lsf/lsf.h"

double lsf_noise_raw_weight(double pscale, double bgnoise) {
	double d = pscale * pscale * bgnoise * bgnoise;
	return d > 0.0 ? 1.0 / d : 0.0;
}

static double wfwhm_formula(const lsf_settings *s, float wfwhm) {
	if (wfwhm <= 0.f)
		return 0.0;
	if (s->fwhm_max <= s->fwhm_min)
		return 1.0;
	double f = wfwhm;
	if (f < s->fwhm_min) f = s->fwhm_min;
	if (f > s->fwhm_max) f = s->fwhm_max;
	double invfwhmax2 = 1.0 / (s->fwhm_max * s->fwhm_max);
	double invdenom = 1.0 / (1.0 / (s->fwhm_min * s->fwhm_min) - invfwhmax2);
	return MAX((1.0 / (f * f) - invfwhmax2) * invdenom, LSF_WEIGHT_FLOOR);
}

static double nbstars_formula(const lsf_settings *s, int nstars) {
	if (s->stars_max <= s->stars_min)
		return 1.0;
	double n = nstars;
	if (n < s->stars_min) n = s->stars_min;
	if (n > s->stars_max) n = s->stars_max;
	double invdenom = 1.0 / (s->stars_max - s->stars_min);
	return MAX((n - s->stars_min) * (n - s->stars_min) * invdenom * invdenom, LSF_WEIGHT_FLOOR);
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
		s->fwhm_min = fmin;
		s->fwhm_max = fmax;
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
		s->stars_min = nmin;
		s->stars_max = nmax;
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
