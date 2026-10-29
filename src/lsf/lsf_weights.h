#ifndef SRC_LSF_LSF_WEIGHTS_H_
#define SRC_LSF_LSF_WEIGHTS_H_

/* LSF frame weights, with constants frozen from the first batch */

#include "lsf/lsf.h"

int lsf_freeze_weights(lsf_settings *s, int nb_frames, int nchan, const double *noise_raw,
		const float *wfwhm, const int *nstars);
double lsf_frame_weight(const lsf_settings *s, int chan, double noise_raw, float wfwhm, int nstars);
double lsf_noise_raw_weight(double pscale, double bgnoise);
gboolean lsf_weight_range_check(const lsf_settings *s, int nb, const float *wfwhm, const int *nstars,
		int *better, int *worse, int *at_floor);
int lsf_switch_to_inverse(lsf_settings *s, int n, const float *hist_wfwhm, const int *hist_nstars,
		const double *hist_weight);

#endif /* SRC_LSF_LSF_WEIGHTS_H_ */
