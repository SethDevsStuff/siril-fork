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

/* LSF accumulator arithmetic (FITS profile v0.4, Sections 10.2, 10.5, 10.6, 10.7).
 * Pure functions, no I/O. */

#include <math.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_cdf.h>

#include "lsf/lsf.h"

/* weighted mean, sum of weights, M2 and sum of squared weights of n samples,
 * two-pass for accuracy */
void lsf_batch_stats(const float *x, const double *w, int n,
		double *mean, double *W, double *M2, double *W2) {
	double sw = 0.0, sw2 = 0.0, swx = 0.0;
	for (int i = 0; i < n; i++) {
		sw += w[i];
		sw2 += w[i] * w[i];
		swx += w[i] * (double) x[i];
	}
	double m = sw > 0.0 ? swx / sw : 0.0;
	double m2 = 0.0;
	for (int i = 0; i < n; i++) {
		double d = (double) x[i] - m;
		m2 += w[i] * d * d;
	}
	*mean = m;
	*W = sw;
	*M2 = m2;
	*W2 = sw2;
}

/* weighted pairwise update (Section 10.2), history a <- a + batch b */
void lsf_merge(float *mean_a, float *w_a, float *m2_a, float *w2_a, WORD *c_a,
		double mean_b, double w_b, double m2_b, double w2_b, int c_b, gboolean *saturated) {
	if (w_b <= 0.0)
		return;
	double wa = *w_a;
	double W = wa + w_b;
	double delta = mean_b - (double) *mean_a;
	*mean_a = (float) ((double) *mean_a + delta * w_b / W);
	*m2_a = (float) ((double) *m2_a + m2_b + delta * delta * wa * w_b / W);
	*w_a = (float) W;
	*w2_a = (float) ((double) *w2_a + w2_b);
	lsf_count_add(c_a, c_b, saturated);
}

/* adds n to a uint16 count that saturates at 65535 and never wraps
 * (Sections 5.2, 5.4) */
void lsf_count_add(WORD *c, int n, gboolean *saturated) {
	if (n <= 0)
		return;
	int v = (int) *c + n;
	if (v > LSF_COUNT_MAX) {
		v = LSF_COUNT_MAX;
		if (saturated)
			*saturated = TRUE;
	}
	*c = (WORD) v;
}

static double phi(double x) {
	return exp(-0.5 * x * x) / sqrt(2.0 * M_PI);
}

/* c = 1/sqrt(v), v the variance of a standard normal clipped to [-kl, kh]
 * (Section 10.6) */
double lsf_clip_correction(double kl, double kh) {
	if (kl <= 0.0 || kh <= 0.0)
		return 1.0;
	double Z = gsl_cdf_ugaussian_P(kh) - gsl_cdf_ugaussian_P(-kl);
	if (Z <= 0.0)
		return 1.0;
	double a = (kl * phi(kl) + kh * phi(kh)) / Z;
	double b = (phi(kl) - phi(kh)) / Z;
	double v = 1.0 - a - b * b;
	if (v <= 0.0)
		return 1.0;
	return 1.0 / sqrt(v);
}

/* predicted variance of a new sample about MEAN (Section 10.7).
 * Returns a negative value where it is undefined. */
double lsf_sigma_pred2(weightingType scheme, double W, double M2, double W2, int N, double w_new) {
	if (W <= 0.0)
		return -1.0;
	if (scheme == NOISE_WEIGHT) {
		if (N < 2 || w_new <= 0.0)
			return -1.0;
		double s02 = M2 / (double) (N - 1);
		return s02 * (1.0 / w_new + 1.0 / W);
	}
	if (W2 <= 0.0)
		return -1.0;
	double neff = W * W / W2;
	if (neff <= LSF_NEFF_MIN)
		return -1.0;
	double s2 = M2 / (W - W2 / W);
	return s2 * (1.0 + W2 / (W * W));
}

/* sigma_hist = margin * max(c * sigma_pred, sigma_frame), sigma_frame alone
 * when sigma_pred is undefined (Section 10.6) */
double lsf_sigma_hist(double sigma_pred2, double c, double sigma_frame, double margin) {
	double s = sigma_frame > 0.0 ? sigma_frame : 0.0;
	if (sigma_pred2 >= 0.0) {
		double sp = c * sqrt(sigma_pred2);
		if (sp > s)
			s = sp;
	}
	return margin * s;
}

/* statistics of Section 10.7 at one pixel and channel: effective sample
 * size, variance of one sample (unit-weight variance for noise weights) and
 * variance of the mean. Each is set to a negative value where undefined. */
void lsf_pixel_stats(weightingType scheme, double W, double M2, double W2, int N,
		double *neff, double *sigma2, double *var_mean) {
	*neff = *sigma2 = *var_mean = -1.0;
	if (W <= 0.0)
		return;
	if (W2 > 0.0)
		*neff = W * W / W2;
	if (scheme == NOISE_WEIGHT) {
		if (N < 2)
			return;
		*sigma2 = M2 / (double) (N - 1);
		*var_mean = *sigma2 / W;
		return;
	}
	if (*neff <= LSF_NEFF_MIN)
		return;
	*sigma2 = M2 / (W - W2 / W);
	*var_mean = *sigma2 * W2 / (W * W);
}

/* value of the Stack (HDU 0) at one pixel: 0 without data, output-normalized
 * when norm, otherwise clamped to [0,1] as Siril's 32-bit stacks are
 * (Sections 3, 4.1) */
float lsf_stack_value(float mean, float w, gboolean norm, double mini, double range) {
	if (!(w > 0.f))
		return 0.f;
	if (norm)
		return (float) (((double) mean - mini) / range);
	return mean < 0.f ? 0.f : (mean > 1.f ? 1.f : mean);
}

/* Adopt (Section 8): the edited Stack disp becomes the new MEAN. A pixel
 * whose Stack value is still the one MEAN gives was not edited and keeps
 * MEAN exactly, which also keeps the values clamped or rounded by the Stack.
 * Edited pixels get MEAN from disp, inverting the output normalization.
 * Pixels without data (W = 0) stay 0. A NaN or infinite edit is not adopted,
 * since every later merge would keep it: those pixels keep MEAN and are
 * counted in nonfinite, which may be NULL. Returns the number of pixels
 * changed. */
size_t lsf_adopt_plane(float *mean, const float *w, const float *disp, size_t n,
		gboolean norm, double mini, double maxi, size_t *nonfinite) {
	double range = maxi - mini;
	size_t changed = 0;
	for (size_t i = 0; i < n; i++) {
		if (!(w[i] > 0.f) || disp[i] == lsf_stack_value(mean[i], w[i], norm, mini, range))
			continue;
		if (!isfinite(disp[i])) {
			if (nonfinite)
				(*nonfinite)++;
			continue;
		}
		mean[i] = norm ? (float) ((double) disp[i] * range + mini) : disp[i];
		changed++;
	}
	return changed;
}

/* min and max for the output normalization: over all channels, where W > 0
 * and the mean is non-zero, as norm_to_0_1_range() does for stacks */
gboolean lsf_outnorm_range(float *const mean[3], float *const W[3], size_t n, int nchan,
		double *mini, double *maxi) {
	float mn = FLT_MAX, mx = -FLT_MAX;
	for (int c = 0; c < nchan; c++) {
		const float *m = mean[c];
		const float *w = W[c];
		for (size_t i = 0; i < n; i++) {
			if (w[i] > 0.f && m[i] != 0.f && isfinite(m[i])) {
				if (m[i] < mn) mn = m[i];
				if (m[i] > mx) mx = m[i];
			}
		}
	}
	if (mn == FLT_MAX || mx <= mn)
		return FALSE;
	*mini = mn;
	*maxi = mx;
	return TRUE;
}

/* Grid expansion (Section 10.5). Footprints are in registration coordinates,
 * where a frame's transform H maps its pixels onto the grid; their y runs
 * opposite to the FITS storage order. box is xmin, xmax, ymin, ymax. */
void lsf_footprint_init(double box[4]) {
	box[0] = box[2] = DBL_MAX;
	box[1] = box[3] = -DBL_MAX;
}

/* adds the corners of an fx x fy frame transformed by h (row-major 3x3) */
void lsf_footprint_add(const double h[9], int fx, int fy, double box[4]) {
	const double cx[4] = { 0.0, fx - 1.0, fx - 1.0, 0.0 }, cy[4] = { 0.0, 0.0, fy - 1.0, fy - 1.0 };
	for (int j = 0; j < 4; j++) {
		double w = h[6] * cx[j] + h[7] * cy[j] + h[8];
		double x = (h[0] * cx[j] + h[1] * cy[j] + h[2]) / w;
		double y = (h[3] * cx[j] + h[4] * cy[j] + h[5]) / w;
		if (x < box[0]) box[0] = x;
		if (x > box[1]) box[1] = x;
		if (y < box[2]) box[2] = y;
		if (y > box[3]) box[3] = y;
	}
}

/* whole pixels of the footprint beyond each edge of an rx x ry grid, in
 * storage order (LSF_PAD_*). A pixel outside the grid is covered only when
 * the footprint reaches its center, so sub-pixel overhangs need nothing. */
gboolean lsf_grid_overhang(const double box[4], int rx, int ry, int need[4]) {
	const double eps = 1e-6;
	double over[4];
	over[LSF_PAD_LEFT] = -box[0];
	over[LSF_PAD_RIGHT] = box[1] - (rx - 1.0);
	over[LSF_PAD_LAST] = -box[2];
	over[LSF_PAD_FIRST] = box[3] - (ry - 1.0);
	gboolean any = FALSE;
	for (int s = 0; s < 4; s++) {
		need[s] = over[s] > 0.0 ? (int) floor(over[s] + eps) : 0;
		any |= need[s] > 0;
	}
	return any;
}

/* each side grows in multiples of step pixels */
void lsf_expand_pads(const int need[4], int step, int pad[4]) {
	for (int s = 0; s < 4; s++)
		pad[s] = need[s] > 0 ? (need[s] + step - 1) / step * step : 0;
}

/* T * h, T being the translation by (dx, dy) in registration coordinates:
 * a transform onto the grid becomes one onto the expanded grid */
void lsf_xfrm_shift(double h[9], double dx, double dy) {
	for (int j = 0; j < 3; j++) {
		h[j] += dx * h[6 + j];
		h[3 + j] += dy * h[6 + j];
	}
}

/* copies an rx x ry plane of elsize-byte pixels, storage order, into dst,
 * which holds the padded plane and is already zeroed */
void lsf_pad_plane(const void *src, void *dst, size_t elsize, int rx, int ry, const int pad[4]) {
	size_t nrx = (size_t) rx + pad[LSF_PAD_LEFT] + pad[LSF_PAD_RIGHT];
	for (int y = 0; y < ry; y++)
		memcpy((char *) dst + (((size_t) y + pad[LSF_PAD_FIRST]) * nrx + pad[LSF_PAD_LEFT]) * elsize,
				(const char *) src + (size_t) y * rx * elsize, (size_t) rx * elsize);
}
