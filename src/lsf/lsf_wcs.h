#ifndef SRC_LSF_LSF_WCS_H_
#define SRC_LSF_LSF_WCS_H_

/* WCS and distortion of a Living Stack: the WCS import of lsf_util
 * (Section 8), the registration of frames through their WCS (Section 10.4)
 * and the distortion model recorded in the ledger (Section 9) */

#include <glib.h>

#include "core/siril.h"
#include "lsf/lsf.h"

struct wcsprm;

int lsf_wcs_read_import(const char *file, const lsf_file *f, fits *out, gchar **why);
double lsf_wcs_sip_shift(struct wcsprm *wcs, int rx, int ry);
gboolean lsf_wcs_same(struct wcsprm *a, struct wcsprm *b, int rx, int ry);
gchar *lsf_wcs_describe(struct wcsprm *wcs, int rx, int ry);
int lsf_wcs_homography(struct wcsprm *img, int rx, int ry, struct wcsprm *ref, int ry_ref,
		transformation_type type, Homography *H);
gchar *lsf_disto_cards(const disto_params *d, fits *frame);
int lsf_register_by_wcs(sequence *seq, int layer, struct wcsprm *ref, int ref_ry, transformation_type type,
		gboolean linear, int max_threads);

#endif /* SRC_LSF_LSF_WCS_H_ */
