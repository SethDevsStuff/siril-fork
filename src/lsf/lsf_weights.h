#ifndef SRC_LSF_LSF_WEIGHTS_H_
#define SRC_LSF_LSF_WEIGHTS_H_

/* LSF frame weights, with constants frozen from the first batch */

#include "lsf/lsf.h"

int lsf_freeze_weights(lsf_settings *s, int nb_frames, int nchan, const double *noise_raw,
		const float *wfwhm, const int *nstars);
double lsf_frame_weight(const lsf_settings *s, int chan, double noise_raw, float wfwhm, int nstars);
double lsf_noise_raw_weight(double pscale, double bgnoise);

#endif /* SRC_LSF_LSF_WEIGHTS_H_ */
