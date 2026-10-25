#ifndef SRC_GUI_LSF_STACKING_H_
#define SRC_GUI_LSF_STACKING_H_

#include <glib.h>

/* the Living Stack modes of the Stacking tab, called by stacking.c */
gboolean stack_lsf_mode(void);
gboolean stack_lsf_append_mode(void);
void stack_set_image_output_name(const gchar *name);
void stack_set_lsf_append(const gchar *path);
void lsf_stacking_start(void);
void lsf_stacking_update_hist_widgets(void);
gboolean lsf_stacking_update_widgets(int stack_method, gboolean can_reframe);
void lsf_stacking_update_append(gboolean listed);
void lsf_stacking_invalidate_checks(void);

#endif
