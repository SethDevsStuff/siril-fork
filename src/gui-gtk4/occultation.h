#ifndef _OCCULTATION_DIALOG_H
#define _OCCULTATION_DIALOG_H

#include <gtk/gtk.h>
GtkWidget *get_occult_dialog();
gpointer occultation_worker(gpointer arg);
gboolean end_occultation_worker(gpointer p);
#endif
