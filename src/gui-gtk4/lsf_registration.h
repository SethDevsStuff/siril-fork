#ifndef SRC_GUI_LSF_REGISTRATION_H_
#define SRC_GUI_LSF_REGISTRATION_H_

#include <glib.h>

#include "registration/registration.h"

/* the Living Stack mode of the Registration tab, called by registration.c */
gboolean lsf_registration_mode(void);
gboolean lsf_registration_drizzle_now(void);
const gchar *lsf_registration_path(void);
gboolean lsf_registration_show_ext_ref(void);
gboolean lsf_registration_update_widgets(void);
int lsf_registration_fill_args(struct registration_args *regargs);

#endif
