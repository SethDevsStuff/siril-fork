#ifndef SRC_GUI_LSF_MANAGER_H_
#define SRC_GUI_LSF_MANAGER_H_

#include <glib.h>

#include "lsf/lsf_commands.h"
#include "lsf/lsf_fits.h"

/* the Living Stack dialog (Tools menu): a file's frames, history, checks,
 * repair, compression, copy, export and maps */
void lsf_manager_open(const gchar *path);

/* what the GUI says of an LSF file, from lsf_read_info */
gchar *lsf_gui_describe(const lsf_info *info);

/* a dialog when an LSF job could not replace its file */
void lsf_gui_report_write(const lsf_job_result *res);

#endif
