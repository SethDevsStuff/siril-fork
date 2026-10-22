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

/* Tests of the LSF accumulator arithmetic (FITS profile v0.4) */

#include <math.h>
#include <string.h>
#include <criterion/criterion.h>

#include "lsf/lsf.h"

static const float samples[] = { 0.12f, 0.15f, 0.11f, 0.14f, 0.13f, 0.18f, 0.10f, 0.16f, 0.12f, 0.17f };
static const double weights[] = { 1.0, 0.5, 2.0, 1.5, 1.0, 0.8, 1.2, 0.9, 1.1, 0.7 };

/* merging two batches equals computing over all samples at once */
Test(lsf_math, merge_equals_one_shot) {
	double m, W, M2, W2;
	lsf_batch_stats(samples, weights, 10, &m, &W, &M2, &W2);

	double ma, wa, m2a, w2a, mb, wb, m2b, w2b;
	lsf_batch_stats(samples, weights, 4, &ma, &wa, &m2a, &w2a);
	lsf_batch_stats(samples + 4, weights + 4, 6, &mb, &wb, &m2b, &w2b);
	float mean = (float) ma, w = (float) wa, m2 = (float) m2a, w2 = (float) w2a;
	WORD count = 4;
	gboolean sat = FALSE;
	lsf_merge(&mean, &w, &m2, &w2, &count, mb, wb, m2b, w2b, 6, &sat);

	cr_assert_float_eq(mean, m, 1e-6);
	cr_assert_float_eq(w, W, 1e-5);
	cr_assert_float_eq(m2, M2, 1e-7);
	cr_assert_float_eq(w2, W2, 1e-5);
	cr_assert_eq(count, 10);
	cr_assert(!sat);
}

Test(lsf_math, merge_into_empty_history) {
	float mean = 0.f, w = 0.f, m2 = 0.f, w2 = 0.f;
	WORD count = 0;
	lsf_merge(&mean, &w, &m2, &w2, &count, 0.25, 3.0, 0.01, 3.0, 3, NULL);
	cr_assert_float_eq(mean, 0.25f, 1e-7);
	cr_assert_float_eq(w, 3.f, 1e-7);
	cr_assert_float_eq(m2, 0.01f, 1e-9);
	cr_assert_eq(count, 3);
}

Test(lsf_math, count_saturates) {
	float mean = 0.5f, w = 1.f, m2 = 0.f, w2 = 1.f;
	WORD count = 65530;
	gboolean sat = FALSE;
	lsf_merge(&mean, &w, &m2, &w2, &count, 0.5, 10.0, 0.0, 10.0, 10, &sat);
	cr_assert_eq(count, 65535);
	cr_assert(sat);
}

/* Section 10.7: c is about 1.014 at +-3 sigma and 1.048 at +-2.5 sigma */
Test(lsf_math, clip_correction) {
	cr_assert_float_eq(lsf_clip_correction(3.0, 3.0), 1.014, 1e-3);
	cr_assert_float_eq(lsf_clip_correction(2.5, 2.5), 1.048, 1e-3);
	cr_assert(lsf_clip_correction(10.0, 10.0) < 1.0000001);
}

/* with N equal weights: s2 = M2 / (N - 1) and sigma_pred2 = s2 (1 + 1/N),
 * for both estimators of Section 10.8 */
Test(lsf_math, sigma_pred_equal_weights) {
	int N = 10;
	double M2 = 0.09, s2 = M2 / (N - 1);
	double rel = lsf_sigma_pred2(NO_WEIGHT, N, M2, N, N, 1.0);
	double prec = lsf_sigma_pred2(NOISE_WEIGHT, N, M2, N, N, 1.0);
	cr_assert_float_eq(rel, s2 * (1.0 + 1.0 / N), 1e-12);
	cr_assert_float_eq(prec, s2 * (1.0 + 1.0 / N), 1e-12);
	/* undefined with a single sample */
	cr_assert(lsf_sigma_pred2(NO_WEIGHT, 1.0, 0.0, 1.0, 1, 1.0) < 0.0);
	cr_assert(lsf_sigma_pred2(NOISE_WEIGHT, 1.0, 0.0, 1.0, 1, 1.0) < 0.0);
}

Test(lsf_math, sigma_hist_floor) {
	/* sigma_frame is a floor */
	cr_assert_float_eq(lsf_sigma_hist(0.0001, 1.0, 0.05, 1.0), 0.05, 1e-12);
	cr_assert_float_eq(lsf_sigma_hist(0.04, 1.0, 0.05, 2.0), 0.4, 1e-12);
	/* undefined sigma_pred: sigma_frame alone */
	cr_assert_float_eq(lsf_sigma_hist(-1.0, 1.014, 0.03, 1.0), 0.03, 1e-12);
}

Test(lsf_math, outnorm_range) {
	float m0[] = { 0.f, 0.2f, 0.5f, 0.9f }, w0[] = { 0.f, 1.f, 1.f, 0.f };
	float *mean[3] = { m0, NULL, NULL }, *W[3] = { w0, NULL, NULL };
	double mini, maxi;
	cr_assert(lsf_outnorm_range(mean, W, 4, 1, &mini, &maxi));
	cr_assert_float_eq(mini, 0.2, 1e-7);
	cr_assert_float_eq(maxi, 0.5, 1e-7);
}

/* the SHA-256 is computed over big-endian float32 values */
Test(lsf_math, hash_of_float_data) {
	float zeros[4] = { 0.f, 0.f, 0.f, 0.f };
	gchar *h = lsf_hash_float_data(zeros, 4);
	/* SHA-256 of 16 zero bytes */
	cr_assert_str_eq(h, "374708fff7719dd5979ec875d56cd2286f6d3cf7ec317a3b25632aab28ec37bb");
	g_free(h);
	float one = 1.f;	// 0x3f800000, big-endian bytes 3f 80 00 00
	gchar *h1 = lsf_hash_float_data(&one, 1);
	guchar be[4] = { 0x3f, 0x80, 0x00, 0x00 };
	gchar *ref = g_compute_checksum_for_data(G_CHECKSUM_SHA256, be, 4);
	cr_assert_str_eq(h1, ref);
	g_free(h1);
	g_free(ref);
}

/* Section 4.1: without output normalization the Stack is MEAN clamped to
 * [0,1], as Siril clamps its 32-bit stacks; 0 where W = 0 */
Test(lsf_math, stack_value) {
	cr_assert_eq(lsf_stack_value(0.25f, 1.f, FALSE, 0.0, 1.0), 0.25f);
	cr_assert_eq(lsf_stack_value(1.57f, 1.f, FALSE, 0.0, 1.0), 1.f);
	cr_assert_eq(lsf_stack_value(-0.01f, 1.f, FALSE, 0.0, 1.0), 0.f);
	cr_assert_eq(lsf_stack_value(0.25f, 0.f, FALSE, 0.0, 1.0), 0.f);
	cr_assert_float_eq(lsf_stack_value(0.6f, 1.f, TRUE, 0.2, 0.8), 0.5f, 1e-7);
}

/* Adopt keeps MEAN exactly where the Stack was not edited (including clamped
 * pixels), takes edited pixels, and leaves pixels without data at 0 */
Test(lsf_math, adopt_unedited_and_edited) {
	float mean[] = { 0.25f, 1.57f, 0.40f, 0.f }, w[] = { 1.f, 1.f, 1.f, 0.f };
	float disp[4];
	for (int i = 0; i < 4; i++)
		disp[i] = lsf_stack_value(mean[i], w[i], FALSE, 0.0, 1.0);
	cr_assert_eq(lsf_adopt_plane(mean, w, disp, 4, FALSE, 0.0, 0.0), 0);
	cr_assert_eq(mean[1], 1.57f);
	disp[2] = 0.30f;	// edited
	disp[3] = 0.50f;	// painted where there is no data
	cr_assert_eq(lsf_adopt_plane(mean, w, disp, 4, FALSE, 0.0, 0.0), 1);
	cr_assert_eq(mean[2], 0.30f);
	cr_assert_eq(mean[3], 0.f);
}

/* with output normalization, an edited pixel is mapped back to MEAN units */
Test(lsf_math, adopt_normalized) {
	double mini = 0.2, maxi = 0.8;
	float mean[] = { 0.2f, 0.5f, 0.8f }, w[] = { 1.f, 1.f, 1.f }, disp[3];
	for (int i = 0; i < 3; i++)
		disp[i] = lsf_stack_value(mean[i], w[i], TRUE, mini, maxi - mini);
	cr_assert_eq(lsf_adopt_plane(mean, w, disp, 3, TRUE, mini, maxi), 0);
	disp[1] = 0.75f;
	cr_assert_eq(lsf_adopt_plane(mean, w, disp, 3, TRUE, mini, maxi), 1);
	cr_assert_float_eq(mean[1], 0.65f, 1e-6);
}

/* Section 10.8 with N equal weights: n_eff = N, s2 = M2 / (N - 1) and the
 * variance of the mean s2 / N, for both estimators */
Test(lsf_math, pixel_stats_equal_weights) {
	int N = 10;
	double M2 = 0.09, s2 = M2 / (N - 1);
	double neff, sig2, vm;
	lsf_pixel_stats(NO_WEIGHT, N, M2, N, N, &neff, &sig2, &vm);
	cr_assert_float_eq(neff, N, 1e-12);
	cr_assert_float_eq(sig2, s2, 1e-12);
	cr_assert_float_eq(vm, s2 / N, 1e-12);
	lsf_pixel_stats(NOISE_WEIGHT, N, M2, N, N, &neff, &sig2, &vm);
	cr_assert_float_eq(neff, N, 1e-12);
	cr_assert_float_eq(sig2, s2, 1e-12);
	cr_assert_float_eq(vm, s2 / N, 1e-12);
	/* one sample: n_eff defined, the variances are not */
	lsf_pixel_stats(NO_WEIGHT, 1.0, 0.0, 1.0, 1, &neff, &sig2, &vm);
	cr_assert_float_eq(neff, 1.0, 1e-12);
	cr_assert(sig2 < 0.0 && vm < 0.0);
	/* no data */
	lsf_pixel_stats(NO_WEIGHT, 0.0, 0.0, 0.0, 0, &neff, &sig2, &vm);
	cr_assert(neff < 0.0 && sig2 < 0.0 && vm < 0.0);
}

/* Section 10.6: a frame shifted +150 columns and +70 storage rows on a
 * 1600x1200 grid has H = translation (150, -70) in registration coordinates
 * (measured on NGC7023 crops). It overhangs the right side and the end of
 * storage, and its shifted transform then fits the expanded grid. */
Test(lsf_math, expansion_shifted_frame) {
	double h[9] = { 1, 0, 150, 0, 1, -70, 0, 0, 1 }, box[4];
	int need[4], pad[4];
	lsf_footprint_init(box);
	lsf_footprint_add(h, 1600, 1200, box);
	cr_assert(lsf_grid_overhang(box, 1600, 1200, need));
	cr_assert_eq(need[LSF_PAD_LEFT], 0);
	cr_assert_eq(need[LSF_PAD_RIGHT], 150);
	cr_assert_eq(need[LSF_PAD_FIRST], 0);
	cr_assert_eq(need[LSF_PAD_LAST], 70);
	lsf_expand_pads(need, 64, pad);
	cr_assert_eq(pad[LSF_PAD_RIGHT], 192);
	cr_assert_eq(pad[LSF_PAD_LAST], 128);
	cr_assert_eq(pad[LSF_PAD_LEFT] + pad[LSF_PAD_FIRST], 0);

	lsf_xfrm_shift(h, pad[LSF_PAD_LEFT], pad[LSF_PAD_LAST]);
	lsf_footprint_init(box);
	lsf_footprint_add(h, 1600, 1200, box);
	cr_assert(!lsf_grid_overhang(box, 1600 + 192, 1200 + 128, need));
}

/* the opposite shift needs columns on the left and rows before row 1 */
Test(lsf_math, expansion_other_side) {
	double h[9] = { 1, 0, -10, 0, 1, 65, 0, 0, 1 }, box[4];
	int need[4], pad[4];
	lsf_footprint_init(box);
	lsf_footprint_add(h, 100, 100, box);
	cr_assert(lsf_grid_overhang(box, 100, 100, need));
	cr_assert_eq(need[LSF_PAD_LEFT], 10);
	cr_assert_eq(need[LSF_PAD_FIRST], 65);
	cr_assert_eq(need[LSF_PAD_RIGHT] + need[LSF_PAD_LAST], 0);
	lsf_expand_pads(need, 64, pad);
	cr_assert_eq(pad[LSF_PAD_LEFT], 64);
	cr_assert_eq(pad[LSF_PAD_FIRST], 128);
}

/* sub-pixel overhangs from registration noise need no expansion */
Test(lsf_math, expansion_subpixel) {
	double h[9] = { 1, 0, -0.4, 0, 1, 0.9, 0, 0, 1 }, box[4];
	int need[4];
	lsf_footprint_init(box);
	lsf_footprint_add(h, 100, 100, box);
	cr_assert(!lsf_grid_overhang(box, 100, 100, need));
	h[2] = -1.0;
	lsf_footprint_init(box);
	lsf_footprint_add(h, 100, 100, box);
	cr_assert(lsf_grid_overhang(box, 100, 100, need));
	cr_assert_eq(need[LSF_PAD_LEFT], 1);
}

/* the shifted transform maps every point to the old image plus the shift */
Test(lsf_math, xfrm_shift) {
	double h[9] = { 0.99, 0.02, 12.5, -0.03, 1.01, -7.25, 1e-6, -2e-6, 1.0 }, s[9];
	memcpy(s, h, sizeof h);
	lsf_xfrm_shift(s, 64.0, 128.0);
	double x = 321.0, y = 123.0;
	double w = h[6] * x + h[7] * y + h[8], ws = s[6] * x + s[7] * y + s[8];
	cr_assert_float_eq((s[0] * x + s[1] * y + s[2]) / ws, (h[0] * x + h[1] * y + h[2]) / w + 64.0, 1e-9);
	cr_assert_float_eq((s[3] * x + s[4] * y + s[5]) / ws, (h[3] * x + h[4] * y + h[5]) / w + 128.0, 1e-9);
}

/* padding keeps every existing pixel and adds zeros (W = 0, no data) */
Test(lsf_math, pad_plane) {
	float src[6] = { 1, 2, 3, 4, 5, 6 };	// 3 x 2, storage order
	int pad[4] = { 1, 2, 1, 0 };
	float dst[6 * 3] = { 0 };
	lsf_pad_plane(src, dst, sizeof(float), 3, 2, pad);
	const float expect[18] = { 0, 0, 0, 0, 0, 0,
		0, 1, 2, 3, 0, 0,
		0, 4, 5, 6, 0, 0 };
	for (int i = 0; i < 18; i++)
		cr_assert_eq(dst[i], expect[i], "pixel %d", i);
	WORD c[2] = { 7, 9 }, cd[2 * 2] = { 0 };
	int p2[4] = { 0, 0, 0, 1 };
	lsf_pad_plane(c, cd, sizeof(WORD), 2, 1, p2);
	cr_assert(cd[0] == 7 && cd[1] == 9 && cd[2] == 0 && cd[3] == 0);
}
