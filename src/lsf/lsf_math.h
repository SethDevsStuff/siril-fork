#ifndef SRC_LSF_LSF_MATH_H_
#define SRC_LSF_LSF_MATH_H_

/* LSF accumulator arithmetic: pure functions, unit-tested */

#include <stddef.h>

#include <glib.h>

#include "core/siril.h"
#include "stacking/stacking.h"

void lsf_batch_stats(const float *x, const double *w, int n,
		double *mean, double *W, double *M2, double *W2);
void lsf_merge(float *mean_a, float *w_a, float *m2_a, float *w2_a, WORD *c_a,
		double mean_b, double w_b, double m2_b, double w2_b, int c_b, gboolean *saturated);
void lsf_count_add(WORD *c, int n, gboolean *saturated);
void lsf_exp_backfill(float *exp, const WORD *count, size_t n, double exptime);
double lsf_clip_correction(double kl, double kh);
double lsf_sigma_pred2(weightingType scheme, double W, double M2, double W2, int N, double w_new);
gboolean lsf_history_testable(double W, double W2);
double lsf_sigma_hist(double sigma_pred2, double c, double sigma_frame, double margin);
gboolean lsf_outnorm_range(float *const mean[3], float *const W[3], size_t n, int nchan,
		double *mini, double *maxi);
void lsf_pixel_stats(weightingType scheme, double W, double M2, double W2, int N,
		double *neff, double *sigma2, double *var_mean);
float lsf_stack_value(float mean, float w, gboolean norm, double mini, double range);
size_t lsf_adopt_plane(float *mean, const float *w, const float *disp, size_t n,
		gboolean norm, double mini, double maxi, size_t *nonfinite);
void lsf_footprint_init(double box[4]);
void lsf_footprint_add(const double h[9], int fx, int fy, double inset, double box[4]);
gboolean lsf_grid_overhang(const double box[4], int rx, int ry, int need[4]);
void lsf_frame_origin(const double h[9], int rx, int ry, int nry, int origin[2]);
int lsf_edge_margin(opencv_interpolation interp);
void lsf_edge_limits(const guint8 *has, int rx, int ry, int margin,
		int *row_lo, int *row_hi, int *col_lo, int *col_hi);
void lsf_expand_pads(const int need[4], int step, int pad[4]);
void lsf_xfrm_shift(double h[9], double dx, double dy);
void lsf_pad_plane(const void *src, void *dst, size_t elsize, int rx, int ry, const int pad[4]);
double lsf_sample_unit(int bitpix);
void lsf_stats_to_unit(imstats *st, double unit);

#endif /* SRC_LSF_LSF_MATH_H_ */
