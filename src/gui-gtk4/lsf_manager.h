#ifndef _GUI_LSF_MANAGER_H_
#define _GUI_LSF_MANAGER_H_

#include <glib.h>
#include "lsf/lsf.h"

/* the Living Stack dialog (Tools menu): a file's frames, history, checks,
 * repair, compression, copy, export and maps (UI plan, B4) */
void lsf_manager_open(const gchar *path);

/* what the GUI says of an LSF file, from lsf_read_info */
gchar *lsf_gui_describe(const lsf_info *info);

#endif
