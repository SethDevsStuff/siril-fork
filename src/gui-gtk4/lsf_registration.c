/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

/* Registration to an existing Living Stack, in the Registration
 * tab: global star alignment against the stacked image of the file,
 * transforms only, with the registration settings stored in it. The external
 * reference of the sequence is set to the file, and lsf_prepare_registration()
 * does the rest, as for the register command. lsf_append then uses these
 * transforms. registration.c calls the lsf_registration_* functions. */

#include <glib/gstdio.h>

#include "core/siril.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/gui_state.h"
#include "gui-gtk4/lsf_manager.h"
#include "gui-gtk4/lsf_registration.h"
#include "gui-gtk4/registration.h"
#include "gui-gtk4/utils.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "lsf/lsf_commands.h"
#include "lsf/lsf_fits.h"

static lsf_info *reg_lsf = NULL;	// the file chosen, read by lsf_read_info, or NULL
static gchar *reg_lsf_problem = NULL;	// why it cannot be used, or NULL
static guint reg_lsf_serial = 0;	// changes each time reg_lsf is read again

static GtkWidget *w(const char *id) {
	return GTK_WIDGET(gtk_builder_get_object(gui.builder, id));
}

gboolean lsf_registration_mode(void) {
	return siril_toggle_get_active(w("reg_lsf_check"));
}

const gchar *lsf_registration_path(void) {
	return gtk_editable_get_text(GTK_EDITABLE(w("reg_lsf_entry")));
}

/* reads the file chosen, and sets the transform and the minimum pairs to
 * those stored in it, which the user may change */
static void load_reg_lsf(void) {
	if (reg_lsf) {
		lsf_info_free(reg_lsf);
		g_free(reg_lsf);
		reg_lsf = NULL;
	}
	reg_lsf_serial++;
	g_clear_pointer(&reg_lsf_problem, g_free);
	gchar *path = g_strstrip(g_strdup(lsf_registration_path()));
	if (!path[0])
		reg_lsf_problem = g_strdup(_("Choose a Living Stack file."));
	else if (!g_file_test(path, G_FILE_TEST_IS_REGULAR))
		reg_lsf_problem = g_strdup(_("This file does not exist."));
	else {
		lsf_info *info = g_new0(lsf_info, 1);
		if (lsf_read_info(path, info)) {
			reg_lsf_problem = g_strdup_printf(_("This is not a Living Stack that can be used: %s."), info->problem);
			lsf_info_free(info);
			g_free(info);
		} else {
			reg_lsf = info;
			if (info->no_append)
				reg_lsf_problem = g_strdup_printf(_("Frames cannot be appended to this file with this version: %s."),
						info->no_append);
			gtk_drop_down_set_selected(GTK_DROP_DOWN(w("comboreg_transfo")), info->file.settings.transform);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(w("spinbut_minpairs")), info->file.settings.min_pairs);
		}
	}
	g_free(path);
	update_reg_interface(TRUE);
}

void on_reg_lsf_check_toggled(GtkCheckButton *button, gpointer user_data) {
	if (gtk_check_button_get_active(button))
		load_reg_lsf();
	else update_reg_interface(TRUE);
}

void on_reg_lsf_entry_changed(GtkEditable *editable, gpointer user_data) {
	load_reg_lsf();
}

void on_reg_lsf_button_clicked(GtkButton *button, gpointer user_data) {
	GtkWindow *control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
	SirilFileChooser *fc = siril_fc_open(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
	siril_fc_set_current_folder_path(fc, com.wd);
	siril_fc_set_select_multiple(fc, FALSE);
	siril_fc_add_filter_pattern(fc, _("Living Stack files (*.lsf.fits)"), "*.lsf.fits;*.LSF.FITS", TRUE);
	if (siril_fc_run(fc) == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_fc_get_filename(fc);
		if (file) {
			GtkEditable *entry = GTK_EDITABLE(w("reg_lsf_entry"));
			gtk_editable_set_text(entry, file);	// loaded by the changed handler
			gtk_editable_set_position(entry, -1);
			g_free(file);
		}
	}
	siril_fc_destroy(fc);
}

/* whether path is a Living Stack, cached on its modification time:
 * update_reg_interface() asks on most widget changes */
static gboolean ref_is_lsf(const gchar *path) {
	static gchar *cached_path = NULL;
	static gint64 cached_mtime = -1, cached_size = -1;
	static gboolean cached = FALSE;
	GStatBuf st;
	if (g_stat(path, &st))
		return FALSE;
	if (g_strcmp0(path, cached_path) || (gint64) st.st_mtime != cached_mtime || (gint64) st.st_size != cached_size) {
		g_free(cached_path);
		cached_path = g_strdup(path);
		cached_mtime = st.st_mtime;
		cached_size = st.st_size;
		cached = fits_file_is_lsf(path);
	}
	return cached;
}

/* A sequence registered against a Living Stack shows it again: sets the
 * file and the mode from the external reference of com.seq. Returns TRUE
 * when it turned the mode on, whose handler comes back to
 * update_reg_interface(), which then returns. */
gboolean lsf_registration_show_ext_ref(void) {
	if (!com.seq.ext_ref || !com.seq.ext_ref_path || !ref_is_lsf(com.seq.ext_ref_path))
		return FALSE;
	GtkEditable *entry = GTK_EDITABLE(w("reg_lsf_entry"));
	if (g_strcmp0(gtk_editable_get_text(entry), com.seq.ext_ref_path))
		gtk_editable_set_text(entry, com.seq.ext_ref_path);	// its handler reads the file
	if (lsf_registration_mode())
		return FALSE;
	siril_toggle_set_active(w("reg_lsf_check"), TRUE);
	return TRUE;
}

/* lsf_check_batch() of the file chosen against com.seq, which reads the
 * reference frame: cached on the file read and the sequence */
static gchar *check_batch_cached(gboolean *fatal) {
	static guint serial = 0;
	static gchar *seqname = NULL, *msg = NULL;
	static int ref = -1, layers = -1;
	static gboolean was_fatal = TRUE;
	int cur_ref = sequence_find_refimage(&com.seq);
	if (serial != reg_lsf_serial || g_strcmp0(seqname, com.seq.seqname) || ref != cur_ref
			|| layers != com.seq.nb_layers) {
		g_free(msg);
		msg = lsf_check_batch(reg_lsf, &com.seq, &was_fatal);
		serial = reg_lsf_serial;
		g_free(seqname);
		seqname = g_strdup(com.seq.seqname);
		ref = cur_ref;
		layers = com.seq.nb_layers;
	}
	*fatal = was_fatal;
	return g_strdup(msg);
}

/* the state of the tab in Living Stack mode, after update_reg_interface()
 * set it for the method: returns whether a registration can start */
gboolean lsf_registration_update_widgets(void) {
	GtkWidget *info = w("reg_lsf_info");
	gboolean lsf = lsf_registration_mode();
	gtk_widget_set_sensitive(w("reg_lsf_box"), lsf);
	gtk_widget_set_visible(info, lsf);
	gtk_widget_set_sensitive(w("comboboxregmethod"), !lsf);
	gtk_widget_set_sensitive(w("reg_2pass"), !lsf);
	if (!lsf) {
		gtk_widget_set_sensitive(w("comboreg_undistort"), TRUE);
		return TRUE;
	}
	/* the external reference of the sequence is the file */
	const gchar *path = lsf_registration_path();
	GtkEditable *ref_entry = GTK_EDITABLE(w("reg_reference_entry"));
	gtk_check_button_set_active(GTK_CHECK_BUTTON(w("reg_reference_checkbutton")), TRUE);
	if (g_strcmp0(gtk_editable_get_text(ref_entry), path))
		gtk_editable_set_text(ref_entry, path);
	gtk_widget_set_sensitive(w("reg_reference_checkbutton"), FALSE);
	gtk_widget_set_sensitive(w("reg_referencefilechooser_box"), FALSE);
	/* no distortion correction, and no output: the append applies */
	GtkDropDown *undistort = GTK_DROP_DOWN(w("comboreg_undistort"));
	if (gtk_drop_down_get_selected(undistort) != DISTO_UNDEF)
		gtk_drop_down_set_selected(undistort, DISTO_UNDEF);
	gtk_widget_set_sensitive(GTK_WIDGET(undistort), FALSE);
	gtk_widget_set_visible(w("output_reg_frame"), FALSE);

	gboolean fatal = reg_lsf == NULL || reg_lsf_problem != NULL;
	gchar *batch = NULL;
	if (reg_lsf && !reg_lsf_problem)
		batch = check_batch_cached(&fatal);
	gchar *summary = reg_lsf ? lsf_gui_describe(reg_lsf) : NULL;
	GString *text = g_string_new(summary ? summary : "");
	if (reg_lsf_problem)
		g_string_append_printf(text, "%s%s", text->len ? "\n" : "", reg_lsf_problem);
	if (batch)
		g_string_append_printf(text, "%s%s", text->len ? "\n" : "", batch);
	if (!fatal)
		g_string_append_printf(text, "%s%s", text->len ? "\n" : "", _("Only the transforms are computed: the append "
					"to this Living Stack applies them. Do not apply them before appending."));
	gtk_label_set_text(GTK_LABEL(info), text->str);
	g_string_free(text, TRUE);
	g_free(summary);
	g_free(batch);
	return !fatal;
}

/* the registration arguments in Living Stack mode: the file as external
 * reference, transforms only, with the star detection stored in the file;
 * the transform and the minimum pairs are those of the widgets, which
 * load_reg_lsf() set from the file. Returns 0, or 1 with the reason logged. */
int lsf_registration_fill_args(struct registration_args *regargs) {
	g_free(regargs->external_ref_path);
	regargs->external_ref_path = g_strdup(lsf_registration_path());
	regargs->use_external_ref = TRUE;
	regargs->undistort = DISTO_UNDEF;
	return lsf_prepare_registration(regargs, TRUE, TRUE, FALSE);
}
