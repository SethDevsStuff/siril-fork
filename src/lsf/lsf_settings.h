#ifndef SRC_LSF_LSF_SETTINGS_H_
#define SRC_LSF_LSF_SETTINGS_H_

/* LSF_SETTINGS: the settings stored in the file, as rows of key, type and value */

#include <glib.h>

#include "core/settings.h"
#include "lsf/lsf.h"
#include "stacking/stacking.h"

void lsf_settings_init(lsf_settings *s);
void lsf_settings_free(lsf_settings *s);
const char *lsf_rejection_name(rejection r);
gboolean lsf_rejection_from_name(const char *name, rejection *r);
const char *lsf_norm_name(normalization n);
gboolean lsf_norm_from_name(const char *name, normalization *n);
const char *lsf_weight_name(weightingType w);
gboolean lsf_weight_from_name(const char *name, weightingType *w);
gboolean lsf_rejection_has_sigma(rejection r);
GPtrArray *lsf_settings_to_rows(const lsf_settings *s, int nchan);
int lsf_settings_from_rows(lsf_settings *s, GPtrArray *rows, int nchan);
void lsf_setting_row_free(gpointer p);

#endif /* SRC_LSF_LSF_SETTINGS_H_ */
