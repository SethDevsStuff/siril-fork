/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
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

#include <stdlib.h>
#include <math.h>

#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/gui_iface.h"
#include "core/processing.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/lsf_manager.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/stacking.h"
#include "gui-gtk4/utils.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "lsf/lsf.h"
#include "registration/mpp.h"
#include "registration/mpp/mpp_config.h"
#include "registration/registration.h"
#include "stacking/stacking.h"
#include "stacking/sum.h"

static struct stacking_args stackparam = { 0 };

static struct filtering_tuple stackfilters[MAX_FILTERS];
/* Values for seq filtering for % or k value*/
static float filter_initvals[] = {90., 3.}; // max %, max k
static float filter_maxvals[] = {100., 5.}; // max %, max k
static float filter_increments[] = {1., 0.1}; // spin button steps for % and k
// update_adjustment passed here as a static (instead of function parameter like in registration)
// in order not to mess up all the calls to update_stack_interface
static int update_adjustment = -1;
static char *filter_tooltip_text[] = {
	N_("Percents of the images of the sequence."),
	N_("Number of standard deviations for rejection of worst images by k-sigma clipping algorithm.")
};

stack_method stacking_methods[] = {
	stack_summing_generic, stack_mean_with_rejection, stack_median, stack_addmax, stack_addmin,
	stack_mpp_handler
};

/* Mpp drizzle scale dropdown. User picks a scale factor; the dispatcher
 * (mpp_stack_apply) then performs all scaling via cv::resize. Entries map
 * one-to-one to the GtkStringList items in siril_stacking.ui (1x / 1.5x /
 * 2x / 3x); sequence-type-specific gating happens at run time. */
static const double mpp_output_scales[] = { 1.0, 1.5, 2.0, 3.0 };

void initialize_stacking_methods() {
	init_stacking_args(&stackparam);
	GtkDropDown *stackcombo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboboxstack_methods"));
	GtkDropDown *rejectioncombo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comborejection"));
	GtkDropDown *weightingcombo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboweighing"));
	GtkSpinButton *low = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
	GtkSpinButton *high = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));
	gtk_drop_down_set_selected(GTK_DROP_DOWN(stackcombo), com.pref.stack.method);
	gtk_drop_down_set_selected(GTK_DROP_DOWN(rejectioncombo), com.pref.stack.rej_method);
	gtk_drop_down_set_selected(GTK_DROP_DOWN(weightingcombo), com.pref.stack.weighting_method);
	switch (gtk_drop_down_get_selected(GTK_DROP_DOWN(rejectioncombo))) {
	case PERCENTILE:
		gtk_spin_button_set_value(low, com.pref.stack.percentile_low);
		gtk_spin_button_set_value(high, com.pref.stack.percentile_high);
		break;
	case LINEARFIT:
		gtk_spin_button_set_value(low, com.pref.stack.linear_low);
		gtk_spin_button_set_value(high, com.pref.stack.linear_high);
		break;
	case SIGMA:
	case MAD:
	case SIGMEDIAN:
	case WINSORIZED:
		gtk_spin_button_set_value(low, com.pref.stack.sigma_low);
		gtk_spin_button_set_value(high, com.pref.stack.sigma_high);
		break;
	case GESDT:
		gtk_spin_button_set_value(low, 0.3);
		gtk_spin_button_set_value(high, 0.05);
		break;
	default:
		return;
	}
}

/* Stacking result: an image, or a Living Stack (LSF) file made by lsf_init
 * from the same widgets, with the LSF options of stack_lsf_grid, or an
 * append to one with lsf_append, which takes its settings from the file */
enum { STACK_RESULT_IMAGE, STACK_RESULT_LSF, STACK_RESULT_LSF_IMAGE, STACK_RESULT_APPEND };

/* stack_lsf_framing items, in order */
static const framing_type lsf_framings[] = { FRAMING_MAX, FRAMING_CURRENT, FRAMING_MIN, FRAMING_COG };

/* The "Save in" entry holds a name pattern in both modes, parsed when the
 * result is written (path_parse). The LSF one starts from the stack-name
 * preference, as a stacking result does, with .lsf.fits (UI plan, C7). */
static gboolean lsf_entry_shown = FALSE;	// the entry holds the LSF pattern
static gchar *other_output = NULL;		// the pattern of the mode not shown
static guint shown_mode = STACK_RESULT_IMAGE;	// the mode the LSF options were last set for

static gchar *default_lsf_pattern() {
	const gchar *pref = com.pref.prepro.stack_default;
	gchar *stem = g_strdup(com.pref.prepro.use_stack_default && pref && pref[0] ? pref : "$seqname$stacked");
	static const char *exts[] = { ".fits", ".fit", ".fts" };
	for (int i = 0; i < (int) G_N_ELEMENTS(exts) && !g_str_has_suffix(stem, ".lsf.fits"); i++)
		if (g_str_has_suffix(stem, exts[i])) {
			stem[strlen(stem) - strlen(exts[i])] = '\0';
			break;
		}
	gchar *pattern = g_str_has_suffix(stem, ".lsf.fits") ? g_strdup(stem) : g_strdup_printf("%s.lsf.fits", stem);
	g_free(stem);
	return pattern;
}

static guint stack_result_mode() {
	static GtkDropDown *mode = NULL;
	if (!mode)
		mode = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_result_mode"));
	return gtk_drop_down_get_selected(mode);
}

/* creating or appending to a Living Stack */
static gboolean stack_lsf_mode() {
	return stack_result_mode() != STACK_RESULT_IMAGE;
}

static gboolean stack_lsf_append_mode() {
	return stack_result_mode() == STACK_RESULT_APPEND;
}

/* the file chosen for an append, read by lsf_read_info, or NULL */
static lsf_info *append_info = NULL;
static gchar *append_summary = NULL;	// what the info label says of it
static void load_append_info(gboolean fill);

/* sets the image mode's stacking name, also while the LSF pattern is shown */
void stack_set_image_output_name(const gchar *name) {
	if (lsf_entry_shown) {
		g_free(other_output);
		other_output = g_strdup(name);
	} else {
		gtk_editable_set_text(GTK_EDITABLE(gtk_builder_get_object(gui.builder, "entryresultfile")), name);
	}
}

void on_stack_result_mode_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)obj;
	(void)pspec;
	GtkEditable *entry = GTK_EDITABLE(gtk_builder_get_object(gui.builder, "entryresultfile"));
	guint mode = stack_result_mode();
	gboolean lsf = mode == STACK_RESULT_LSF || mode == STACK_RESULT_LSF_IMAGE;
	/* a new file keeps rejection maps by default (UI plan, C7); an append
	 * adds them to a file without only on request */
	if ((mode == STACK_RESULT_APPEND) != (shown_mode == STACK_RESULT_APPEND))
		siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps")),
				mode != STACK_RESULT_APPEND);
	shown_mode = mode;
	if (mode == STACK_RESULT_APPEND) {
		/* the entry is hidden, the file is chosen in stack_lsf_append_entry. The
		 * widgets may have been changed in another mode: show the file again */
		load_append_info(TRUE);
		return;
	} else if (lsf != lsf_entry_shown) {
		/* swap the patterns of the two modes */
		gchar *shown = g_strdup(gtk_editable_get_text(entry));
		if (other_output)
			gtk_editable_set_text(entry, other_output);
		else {
			gchar *pattern = default_lsf_pattern();
			gtk_editable_set_text(entry, pattern);
			g_free(pattern);
		}
		g_free(other_output);
		other_output = shown;
		lsf_entry_shown = lsf;
		gtk_widget_set_tooltip_text(GTK_WIDGET(entry), lsf ? _("Name of the Living Stack file. Like the name of "
					"a stacking result, it may use $seqname$ and header keys such as $FILTER:%s$, taken from "
					"the reference frame. The .lsf.fits extension is added when it is missing.") : NULL);
	}
	update_stack_interface(TRUE);
}

/* the LSF file name from the pattern, parsed like the name of a stacking
 * result but with the header of the reference frame, since the Stack does
 * not exist yet, and with the .lsf.fits extension. NULL if it is empty. */
static gchar *parse_lsf_name(const gchar *pattern) {
	gchar *expression = g_strstrip(g_strdup(pattern));
	gchar *name = NULL;
	fits ref = { 0 };
	int status = PATHPARSE_ERR_OK;
	if (g_utf8_strchr(expression, -1, '$') &&
			!seq_read_frame_metadata(&com.seq, sequence_find_refimage(&com.seq), &ref))
		name = update_header_and_parse(&ref, expression, PATHPARSE_MODE_WRITE_NOFAIL, TRUE, &status);
	else name = g_strdup(expression);
	clearfits(&ref);
	g_free(expression);
	if (!name || !name[0]) {
		g_free(name);
		return NULL;
	}
	if (g_str_has_suffix(name, ".lsf.fits"))
		return name;
	/* a FITS extension, then separators left by $seqname$ or a missing key */
	static const char *exts[] = { ".fits", ".fit", ".fts" };
	for (int i = 0; i < (int) G_N_ELEMENTS(exts); i++)
		if (g_str_has_suffix(name, exts[i])) {
			name[strlen(name) - strlen(exts[i])] = '\0';
			break;
		}
	size_t len = strlen(name);
	while (len > 1 && (name[len - 1] == '_' || name[len - 1] == '-' || name[len - 1] == '.'))
		name[--len] = '\0';
	gchar *full = g_strdup_printf("%s.lsf.fits", name);
	g_free(name);
	return full;
}

/* the history thresholds lsf_init takes from the rejection (Section 10.6) */
static void lsf_hist_from_rejection(double *low, double *high) {
	GtkDropDown *rejec = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comborejection"));
	GtkSpinButton *siglow = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
	GtkSpinButton *sighigh = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));
	gboolean sigma = lsf_rejection_has_sigma(gtk_drop_down_get_selected(rejec));
	double l = gtk_spin_button_get_value(siglow), h = gtk_spin_button_get_value(sighigh);
	*low = sigma && l > 0.0 ? l : 3.0;
	*high = sigma && h > 0.0 ? h : 3.0;
}

static void update_lsf_hist_widgets() {
	GtkWidget *same = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_hist_same"));
	GtkSpinButton *low = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_low"));
	GtkSpinButton *high = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_high"));
	/* an append starts from the thresholds stored in the file */
	gboolean from_rejection = !stack_lsf_append_mode() && siril_toggle_get_active(same);
	gtk_widget_set_sensitive(GTK_WIDGET(low), !from_rejection);
	gtk_widget_set_sensitive(GTK_WIDGET(high), !from_rejection);
	if (from_rejection) {
		double l, h;
		lsf_hist_from_rejection(&l, &h);
		gtk_spin_button_set_value(low, l);
		gtk_spin_button_set_value(high, h);
	}
}

void on_stack_lsf_hist_same_toggled(GtkCheckButton *button, gpointer user_data) {
	update_lsf_hist_widgets();
}

/* shows the new Stack, which lsf_load_result put in gfit, as stacking does */
static void end_lsf_stacking(const lsf_job_result *res) {
	set_cursor_waiting(FALSE);
	lsf_gui_report_write(res);
	if (res->retval) {
		siril_log_error(_("The Living Stack was not created, please check the log to fix your issue.\n"));
		return;
	}
	gui_iface.on_stack_complete();
}

/* the widgets of the method expander show the settings of the file: the
 * facts (normalization, weighting, output normalization) cannot change and
 * the rejection cannot be overridden yet (UI plan, C2). Filling them goes
 * through their handlers, which would store the file's values as the user's
 * stacking preferences, so these are kept. */
static void fill_widgets_from_lsf(const lsf_settings *s) {
	struct stack_config saved = com.pref.stack;
	gtk_drop_down_set_selected(GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comborejection")), s->rej);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_siglow_button")), s->sig[0]);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_sighigh_button")), s->sig[1]);
	gtk_drop_down_set_selected(GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combonormalize")), s->norm);
	siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkfastnorm")), s->lite_norm);
	siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_RGBequal")), s->rgb_equal);
	gtk_drop_down_set_selected(GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboweighing")), s->weighting);
	siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_normalise_to_max")), s->outnorm);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_low")), s->hist_low);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_high")), s->hist_high);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_margin")), s->hist_margin);
	siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps")), FALSE);
	com.pref.stack = saved;
}

/* the name suggested for the new file of an append: the file's, with the
 * generation it will have, replacing the one an earlier suggestion added */
static gchar *suggest_append_out(const lsf_info *info) {
	const char *path = info->file.path;
	gchar *stem = g_str_has_suffix(path, ".lsf.fits") ? g_strndup(path, strlen(path) - strlen(".lsf.fits")) : g_strdup(path);
	gchar *gen = g_strrstr(stem, "_gen");
	if (gen && gen[4] && strspn(gen + 4, "0123456789") == strlen(gen + 4))
		*gen = '\0';
	gchar *name = g_strdup_printf("%s_gen%d.lsf.fits", stem, info->file.gen + 1);
	g_free(stem);
	return name;
}

/* the file an append writes instead of overwriting the chosen one, or NULL */
static gchar *append_out_path() {
	if (!siril_toggle_get_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_newfile"))))
		return NULL;
	return parse_lsf_name(gtk_editable_get_text(GTK_EDITABLE(gtk_builder_get_object(gui.builder, "stack_lsf_append_out"))));
}

/* reads the file of stack_lsf_append_entry and shows it; with fill, also
 * sets the widgets of the method expander to its settings and suggests a
 * name for a new file */
static void load_append_info(gboolean fill) {
	GtkEditable *entry = GTK_EDITABLE(gtk_builder_get_object(gui.builder, "stack_lsf_append_entry"));
	if (append_info) {
		lsf_info_free(append_info);
		g_free(append_info);
		append_info = NULL;
	}
	g_free(append_summary);
	append_summary = NULL;
	gchar *path = g_strstrip(g_strdup(gtk_editable_get_text(entry)));
	if (!path[0])
		append_summary = g_strdup(_("Choose a Living Stack file."));
	else if (!g_file_test(path, G_FILE_TEST_IS_REGULAR))
		append_summary = g_strdup(_("This file does not exist."));
	else {
		lsf_info *info = g_new0(lsf_info, 1);
		if (lsf_read_info(path, info)) {
			append_summary = g_strdup_printf(_("This file cannot be appended to: %s."), info->problem);
			lsf_info_free(info);
			g_free(info);
		} else {
			append_info = info;
			append_summary = lsf_gui_describe(info);
			if (fill) {
				fill_widgets_from_lsf(&info->file.settings);
				gchar *out = suggest_append_out(info);
				gtk_editable_set_text(GTK_EDITABLE(gtk_builder_get_object(gui.builder, "stack_lsf_append_out")), out);
				g_free(out);
			}
		}
	}
	g_free(path);
	update_stack_interface(TRUE);
}

void on_stack_lsf_append_entry_changed(GtkEditable *editable, gpointer user_data) {
	load_append_info(TRUE);
}

void on_stack_lsf_append_option_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_stack_lsf_append_button_clicked(GtkButton *button, gpointer user_data) {
	GtkWindow *control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
	GtkEditable *entry = GTK_EDITABLE(gtk_builder_get_object(gui.builder, "stack_lsf_append_entry"));
	SirilFileChooser *fc = siril_fc_open(control_window, GTK_FILE_CHOOSER_ACTION_OPEN);
	siril_fc_set_current_folder_path(fc, com.wd);
	siril_fc_set_select_multiple(fc, FALSE);
	siril_fc_add_filter_pattern(fc, _("Living Stack files (*.lsf.fits)"), "*.lsf.fits;*.LSF.FITS", TRUE);
	if (siril_fc_run(fc) == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_fc_get_filename(fc);
		if (file) {
			gtk_editable_set_text(entry, file);	// loaded by the changed handler
			gtk_editable_set_position(entry, -1);
			g_free(file);
		}
	}
	siril_fc_destroy(fc);
}

/* the new file chosen with the file chooser, which already asked before
 * replacing it */
static gchar *append_out_chosen = NULL;

void on_stack_lsf_append_newfile_toggled(GtkCheckButton *button, gpointer user_data) {
	gboolean on = gtk_check_button_get_active(button);
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_out")), on);
	gtk_widget_set_sensitive(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_out_button")), on);
	update_stack_interface(TRUE);
}

void on_stack_lsf_append_out_changed(GtkEditable *editable, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_stack_lsf_append_out_button_clicked(GtkButton *button, gpointer user_data) {
	GtkWindow *control_window = GTK_WINDOW(GTK_APPLICATION_WINDOW(gtk_builder_get_object(gui.builder, "control_window")));
	GtkEditable *entry = GTK_EDITABLE(gtk_builder_get_object(gui.builder, "stack_lsf_append_out"));
	SirilFileChooser *fc = siril_fc_save(control_window, GTK_FILE_CHOOSER_ACTION_SAVE);
	siril_fc_set_current_folder_path(fc, com.wd);
	siril_fc_add_filter_pattern(fc, _("Living Stack files (*.lsf.fits)"), "*.lsf.fits;*.LSF.FITS", TRUE);
	gchar *suggested = g_path_get_basename(gtk_editable_get_text(entry));
	if (suggested[0] && strcmp(suggested, "."))
		siril_fc_set_current_name(fc, suggested);
	g_free(suggested);
	if (siril_fc_run(fc) == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_fc_get_filename(fc);
		if (file) {
			g_free(append_out_chosen);
			append_out_chosen = parse_lsf_name(file);
			gtk_editable_set_text(entry, file);
			gtk_editable_set_position(entry, -1);
			g_free(file);
		}
	}
	siril_fc_destroy(fc);
}

/* after a registration of the sequence against a Living Stack (UI plan, B3):
 * the next stacking appends to it */
void stack_set_lsf_append(const gchar *path) {
	GtkEditable *entry = GTK_EDITABLE(gtk_builder_get_object(gui.builder, "stack_lsf_append_entry"));
	GtkDropDown *mode = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_result_mode"));
	gboolean same_file = !g_strcmp0(gtk_editable_get_text(entry), path);
	if (!same_file) {
		gtk_editable_set_text(entry, path);	// loaded by the changed handler
		gtk_editable_set_position(entry, -1);
	}
	if (stack_result_mode() != STACK_RESULT_APPEND)
		gtk_drop_down_set_selected(mode, STACK_RESULT_APPEND);	// loaded by the mode handler
	else if (same_file)
		load_append_info(TRUE);	// the file's binding note changed
}

/* the problems update_stack_interface shows and start_lsf_append checks:
 * the file, the batch and the pre-flight. Returns a message or NULL. */
static gchar *append_problems(gboolean *fatal) {
	*fatal = TRUE;
	if (!append_info)
		return NULL;
	if (append_info->no_append)
		return g_strdup_printf(_("This version cannot append to this file: %s."), append_info->no_append);
	gboolean batch_fatal = FALSE, flight_fatal = FALSE;
	gchar *batch = lsf_check_batch(append_info, &com.seq, &batch_fatal);
	if (batch_fatal)
		return batch;
	GtkWidget *noreg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_noreg"));
	GtkWidget *rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps"));
	gchar *out = append_out_path();
	gchar *flight = lsf_append_preflight(append_info, com.wd, out, stackparam.nb_images_to_stack, siril_toggle_get_active(noreg),
			gtk_widget_get_visible(rejmaps) && siril_toggle_get_active(rejmaps), &flight_fatal);
	g_free(out);
	*fatal = flight_fatal;
	if (flight_fatal) {
		g_free(batch);
		return flight;
	}
	gchar *msg = batch && flight ? g_strdup_printf("%s %s", batch, flight) : g_strdup(batch ? batch : flight);
	g_free(batch);
	g_free(flight);
	return msg;
}

/* a registration of the sequence against the file, made in the Registration
 * tab or with register -extref= (UI plan, B3): what the append does with it.
 * Returns a note for the info label or NULL; *reusable when it is used. */
static gchar *append_registration_note(gboolean *reusable) {
	*reusable = FALSE;
	if (!append_info)
		return NULL;
	int layer;
	gchar *date = NULL, *why = NULL, *msg = NULL;
	*reusable = lsf_binding_check(&com.seq, &append_info->file, &layer, &date, &why);
	if (*reusable)
		msg = g_strdup_printf(_("The sequence was registered against this Living Stack (%s): the append uses "
					"that registration."), date ? date : "");
	else if (why)
		msg = g_strdup_printf(_("Earlier registration not used, as %s: the append registers the frames again."), why);
	g_free(date);
	g_free(why);
	return msg;
}

/* shows the Stack after the append and the file as it is now */
static void end_lsf_append(const lsf_job_result *res) {
	set_cursor_waiting(FALSE);
	lsf_gui_report_write(res);
	if (res->retval)
		siril_log_error(_("The frames were not appended, please check the log to fix your issue.\n"));
	else gui_iface.on_stack_complete();
	if (!res->retval && res->out) {
		/* as after Save As, the next append goes to the new file */
		siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_newfile")), FALSE);
		g_clear_pointer(&append_out_chosen, g_free);
		gchar *path = g_canonicalize_filename(res->out, NULL);
		stack_set_lsf_append(path);
		g_free(path);
		return;
	}
	load_append_info(TRUE);
}

/* the frames the stacking filters chose (UI plan, A4), for the LSF job, or
 * FALSE when there are none */
static gboolean set_lsf_frames(struct lsf_cmd_args *a) {
	if (!stackparam.image_indices || stackparam.nb_images_to_stack < 1)
		return FALSE;
	a->nb_images = stackparam.nb_images_to_stack;
	a->image_indices = g_memdup2(stackparam.image_indices, a->nb_images * sizeof(int));
	return TRUE;
}

/* the job and the command it is logged as (a command has no stacking filters) */
static void log_and_start_lsf(struct lsf_cmd_args *a) {
	gchar *cmd = lsf_cmd_args_to_command(a);
	if (cmd) {
		gchar **lines = g_strsplit(cmd, "\n", -1);
		for (int i = 0; lines[i]; i++)
			siril_log_message(_("Living Stack, equivalent command: %s\n"), lines[i]);
		g_strfreev(lines);
		g_free(cmd);
	}
	if (a->image_indices && a->nb_images != com.seq.selnum)
		siril_log_message(_("Living Stack: %d frames chosen by the stacking filters; the command uses the "
					"frames selected in the sequence instead\n"), a->nb_images);
	set_cursor_waiting(TRUE);
	if (lsf_start(a)) {
		PRINT_ANOTHER_THREAD_RUNNING;
		set_cursor_waiting(FALSE);
	}
}

/* appends the selected frames to the chosen Living Stack, with lsf_append */
static void start_lsf_append() {
	GtkDropDown *framing_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_lsf_append_framing"));
	GtkDropDown *comp_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_lsf_append_compression"));
	GtkSpinButton *step = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_expand_step"));
	GtkSpinButton *hist_low = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_low"));
	GtkSpinButton *hist_high = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_high"));
	GtkSpinButton *hist_margin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_margin"));
	GtkWidget *force = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_expand_force"));
	GtkWidget *noreg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_noreg"));
	GtkWidget *rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps"));
	GtkWidget *newfile = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_newfile"));
	GtkWidget *export_image = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_export"));
	GtkDropDown *rejec_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comborejection"));
	GtkSpinButton *siglow = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
	GtkSpinButton *sighigh = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}
	/* the user's choices, which may override the defaults stored in the file */
	rejection rej = gtk_drop_down_get_selected(rejec_combo);
	double sig_low = gtk_spin_button_get_value(siglow), sig_high = gtk_spin_button_get_value(sighigh);
	double low = gtk_spin_button_get_value(hist_low), high = gtk_spin_button_get_value(hist_high),
	       margin = gtk_spin_button_get_value(hist_margin);
	/* the file may have changed since it was chosen: read it again, but
	 * keep the widgets as the user set them */
	load_append_info(FALSE);
	if (!append_info) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No Living Stack"), append_summary);
		return;
	}
	/* a new file, or the chosen one overwritten */
	gchar *out = append_out_path();
	if (siril_toggle_get_active(newfile) && !out) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No file name"), _("Enter the name of the new Living Stack file, "
					"or uncheck \"Save to a new file\" to overwrite the chosen one."));
		return;
	}
	if (out && lsf_same_file(out, append_info->file.path))
		g_clear_pointer(&out, g_free);	// the chosen file itself
	gboolean overwrite = FALSE;
	if (out && g_file_test(out, G_FILE_TEST_EXISTS)) {
		gchar *msg = g_strdup_printf(_("%s already exists. It is replaced by the chosen Living Stack with the frames "
					"appended."), out);
		gboolean go = !g_strcmp0(out, append_out_chosen) ||
			siril_confirm_dialog(_("Replace the file?"), msg, _("Replace"));
		g_free(msg);
		if (!go) {
			g_free(out);
			return;
		}
		overwrite = TRUE;
	}
	gboolean fatal;
	gchar *problems = append_problems(&fatal);
	if (problems) {
		gboolean go = !fatal && siril_confirm_dialog(_("Append anyway?"), problems, _("Append"));
		if (fatal)
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Cannot append"), problems);
		g_free(problems);
		if (!go) {
			g_free(out);
			return;
		}
	}

	/* overrides of the rejection (UI plan, C2) and of the history test (B2)
	 * replace the settings stored in the file: compared at the precision
	 * the spin buttons show, and confirmed */
	const lsf_settings *s = &append_info->file.settings;
	gboolean rej_set = rej != s->rej || (rej != NO_REJEC &&
			(fabs(sig_low - s->sig[0]) > 0.005 || fabs(sig_high - s->sig[1]) > 0.005));
	gboolean low_set = fabs(low - s->hist_low) > 0.005, high_set = fabs(high - s->hist_high) > 0.005,
		 margin_set = fabs(margin - s->hist_margin) > 0.005;
	if (rej_set || low_set || high_set || margin_set) {
		GString *changes = g_string_new(NULL);
		if (rej_set) {
			gchar *from = s->rej == NO_REJEC ? g_strdup(lsf_rejection_name(s->rej)) :
				g_strdup_printf("%s %g %g", lsf_rejection_name(s->rej), s->sig[0], s->sig[1]);
			gchar *to = rej == NO_REJEC ? g_strdup(lsf_rejection_name(rej)) :
				g_strdup_printf("%s %g %g", lsf_rejection_name(rej), sig_low, sig_high);
			g_string_append_printf(changes, _("Rejection: from %s to %s.\n"), from, to);
			g_free(from);
			g_free(to);
		}
		if (low_set || high_set || margin_set)
			g_string_append_printf(changes, _("History test: from low %g, high %g, margin %g to low %g, high %g, "
						"margin %g.\n"), s->hist_low, s->hist_high, s->hist_margin, low_set ? low : s->hist_low,
					high_set ? high : s->hist_high, margin_set ? margin : s->hist_margin);
		g_string_append(changes, _("These settings are stored in the Living Stack: the change overwrites them "
					"for this append and the appends that follow. The frames already in the stack are not "
					"affected."));
		gboolean go = siril_confirm_dialog(_("Change the settings of the Living Stack?"), changes->str, _("Append"));
		g_string_free(changes, TRUE);
		if (!go) {
			g_free(out);
			return;
		}
	}

	/* lsf_append reads the sequence from its file: saved with the state
	 * shown, the frames are those chosen by the stacking filters */
	writeseqfile(&com.seq);

	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_APPEND);
	if (!a) {
		g_free(out);
		return;
	}
	a->out = out;	// freed with a
	a->overwrite = overwrite;
	if (!set_lsf_frames(a)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No frames"), _("The stacking filters leave no frame to append."));
		lsf_cmd_args_free(a);
		return;
	}
	a->path = g_strdup(append_info->file.path);
	if (rej_set) {
		a->settings.rej = rej;
		a->settings.sig[0] = rej == NO_REJEC ? 0.0 : sig_low;
		a->settings.sig[1] = rej == NO_REJEC ? 0.0 : sig_high;
		a->rej_set = TRUE;
	}
	a->seqname = g_strdup(com.seq.seqname);
	a->noreg = siril_toggle_get_active(noreg);
	/* shown when a registration against the file can be reused (UI plan, B3) */
	GtkWidget *force_reg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_force_reg"));
	a->force_reg = !a->noreg && gtk_widget_get_visible(force_reg) && siril_toggle_get_active(force_reg);
	a->framing = gtk_drop_down_get_selected(framing_combo) == 1 ? FRAMING_CURRENT : FRAMING_MAX;
	a->framing_set = a->framing != FRAMING_MAX;
	a->expand_step = gtk_spin_button_get_value_as_int(step);
	a->expand_force = siril_toggle_get_active(force);
	a->hist_low = low;
	a->hist_low_set = low_set;
	a->hist_high = high;
	a->hist_high_set = high_set;
	a->hist_margin = margin;
	a->hist_margin_set = margin_set;
	a->rejmaps = gtk_widget_get_visible(rejmaps) && siril_toggle_get_active(rejmaps);
	/* stack_lsf_append_compression: keep, none, gzip2, gzip1 */
	guint comp = gtk_drop_down_get_selected(comp_combo);
	if (comp >= 1 && comp <= 3) {
		a->compression = comp == 1 ? LSF_COMP_NONE : comp == 2 ? LSF_COMP_GZIP2 : LSF_COMP_GZIP1;
		a->compression_set = TRUE;
	}
	/* the plain image of the updated stack, named like the file written (UI plan, C10) */
	const char *target = a->out ? a->out : a->path;
	if (siril_toggle_get_active(export_image) && g_str_has_suffix(target, ".lsf.fits"))
		a->export_out = g_strndup(target, strlen(target) - strlen(".lsf.fits"));
	a->load_result = TRUE;
	a->end_func = end_lsf_append;
	log_and_start_lsf(a);
}

/* creates a Living Stack from the selected frames, with lsf_init */
static void start_lsf_stacking() {
	GtkDropDown *rejec_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comborejection"));
	GtkDropDown *norm_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combonormalize"));
	GtkDropDown *weighing_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboweighing"));
	GtkDropDown *mode_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_result_mode"));
	GtkDropDown *framing_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_lsf_framing"));
	GtkDropDown *comp_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_lsf_compression"));
	GtkSpinButton *siglow = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
	GtkSpinButton *sighigh = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));
	GtkSpinButton *hist_low = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_low"));
	GtkSpinButton *hist_high = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_high"));
	GtkSpinButton *hist_margin = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_margin"));
	GtkWidget *hist_same = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_hist_same"));
	GtkWidget *fast_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkfastnorm"));
	GtkWidget *rgb_equal = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_RGBequal"));
	GtkWidget *output_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_normalise_to_max"));
	GtkWidget *rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps"));
	GtkWidget *overwrite = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkbutoverwrite"));
	GtkEditable *output_file = GTK_EDITABLE(gtk_builder_get_object(gui.builder, "entryresultfile"));

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}
	gchar *path = parse_lsf_name(gtk_editable_get_text(output_file));
	if (!path) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No file name"), _("Enter the name of the Living Stack file to create."));
		return;
	}
	siril_log_message(_("Living Stack file: %s\n"), path);

	normalization norm = gtk_drop_down_get_selected(norm_combo);
	weightingType weighting = gtk_widget_get_visible(GTK_WIDGET(weighing_combo)) ?
		gtk_drop_down_get_selected(weighing_combo) : NO_WEIGHT;
	if (weighting == NBSTACK_WEIGHT) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Weighting not supported"),
				_("A Living Stack cannot be weighted by the number of images. Choose another weighting."));
		g_free(path);
		return;
	}
	if (weighting == NOISE_WEIGHT && norm == NO_NORM) {
		siril_log_error(_("Weighting by noise is allowed only if normalization has been activated, ignoring weights.\n"));
		weighting = NO_WEIGHT;
	}

	/* memory and disk first (UI plan, A6) */
	gboolean fatal = FALSE;
	gchar *flight = lsf_init_preflight(&com.seq, com.wd, path, stackparam.nb_images_to_stack,
			siril_toggle_get_active(rejmaps), &fatal);
	if (flight) {
		gboolean go = !fatal && siril_confirm_dialog(_("Create anyway?"), flight, _("Create"));
		if (fatal)
			siril_message_dialog(GTK_MESSAGE_ERROR, _("Cannot create the Living Stack"), flight);
		g_free(flight);
		if (!go) {
			g_free(path);
			return;
		}
	}

	if (stackparam.nb_images_to_stack < LSF_LOW_STACK_WARN) {
		gchar *msg = g_strdup_printf(_("Only %d frames are selected, fewer than %d. The scatter of each pixel estimated "
					"from so few frames is unreliable, and any rejection decided on it stays in the stack."),
				stackparam.nb_images_to_stack, LSF_LOW_STACK_WARN);
		gboolean go = siril_confirm_dialog(_("Few frames for a Living Stack"), msg, _("Create anyway"));
		g_free(msg);
		if (!go) {
			g_free(path);
			return;
		}
	}
	gboolean replace = siril_toggle_get_active(overwrite);
	if (g_file_test(path, G_FILE_TEST_EXISTS)) {
		if (!replace || !fits_file_is_lsf(path)) {
			gchar *msg = replace ? g_strdup_printf(_("%s exists and is not a Living Stack file. Choose another name."), path) :
				g_strdup_printf(_("%s already exists. Choose another name, or check Overwrite to replace it."), path);
			siril_message_dialog(GTK_MESSAGE_ERROR, _("File exists"), msg);
			g_free(msg);
			g_free(path);
			return;
		}
		gchar *msg = g_strdup_printf(_("%s is a Living Stack already. Replacing it loses every frame accumulated "
					"in it, and there is no undo."), path);
		gboolean go = siril_confirm_dialog(_("Replace the Living Stack?"), msg, _("Replace"));
		g_free(msg);
		if (!go) {
			g_free(path);
			return;
		}
	}

	/* lsf_init reads the sequence from its file: saved with the state shown,
	 * the frames are those chosen by the stacking filters */
	writeseqfile(&com.seq);

	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_INIT);
	if (!a) {
		g_free(path);
		return;
	}
	if (!set_lsf_frames(a)) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("No frames"), _("The stacking filters leave no frame to stack."));
		g_free(path);
		lsf_cmd_args_free(a);
		return;
	}
	lsf_settings *s = &a->settings;
	a->seqname = g_strdup(com.seq.seqname);
	a->path = path;
	s->rej = gtk_drop_down_get_selected(rejec_combo);
	if (s->rej != NO_REJEC) {
		s->sig[0] = gtk_spin_button_get_value(siglow);
		s->sig[1] = gtk_spin_button_get_value(sighigh);
	}
	s->norm = norm;
	s->lite_norm = norm != NO_NORM && siril_toggle_get_active(fast_norm);
	s->rgb_equal = norm != NO_NORM && com.seq.nb_layers == 3 && siril_toggle_get_active(rgb_equal);
	s->weighting = weighting;
	s->outnorm = siril_toggle_get_active(output_norm);
	guint fr = gtk_drop_down_get_selected(framing_combo);
	a->framing = fr < G_N_ELEMENTS(lsf_framings) ? lsf_framings[fr] : FRAMING_MAX;
	a->framing_set = a->framing != FRAMING_MAX;
	if (!siril_toggle_get_active(hist_same)) {
		a->hist_low = gtk_spin_button_get_value(hist_low);
		a->hist_high = gtk_spin_button_get_value(hist_high);
		a->hist_low_set = a->hist_high_set = TRUE;
	}
	if (gtk_spin_button_get_value(hist_margin) != 1.0) {
		a->hist_margin = gtk_spin_button_get_value(hist_margin);
		a->hist_margin_set = TRUE;
	}
	a->rejmaps = siril_toggle_get_active(rejmaps);
	guint comp = gtk_drop_down_get_selected(comp_combo);
	if (comp == 1 || comp == 2) {
		a->compression = comp == 1 ? LSF_COMP_GZIP2 : LSF_COMP_GZIP1;
		a->compression_set = TRUE;
	}
	a->overwrite = replace;
	if (gtk_drop_down_get_selected(mode_combo) == STACK_RESULT_LSF_IMAGE)
		a->export_out = g_strndup(path, strlen(path) - strlen(".lsf.fits"));
	a->load_result = TRUE;
	a->end_func = end_lsf_stacking;
	log_and_start_lsf(a);
}

/* starts a summing operation using data stored in the stackparam structure
 * function is not reentrant but can be called again after it has returned and the thread is running */
static void start_stacking() {
	gchar *error = NULL;
	static GtkDropDown *method_combo = NULL, *rejec_combo = NULL, *norm_combo = NULL, *weighing_combo;
	static GtkEntry *output_file = NULL;
	static GtkCheckButton *overwrite = NULL, *force_norm = NULL, *max_framing = NULL,
					*fast_norm = NULL, *rejmaps = NULL, *merge_rejmaps = NULL, 
					*upscale_at_stacking = NULL, *overlap_norm = NULL, *force32b = NULL;
	static GtkSpinButton *sigSpin[2] = {NULL, NULL}, *feather_dist = NULL;
	static GtkWidget *norm_to_max = NULL, *RGB_equal = NULL, *blend_frame = NULL;
	static GtkDropDown *mpp_scale_combo = NULL, *mpp_engine_combo = NULL;
	static GtkSpinButton *mpp_stack_percent = NULL, *mpp_stack_frames = NULL,
	                     *mpp_bg_fraction = NULL, *mpp_bg_blend = NULL;
	static GtkCheckButton *mpp_skip_failed = NULL;

	if (method_combo == NULL) {
		method_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboboxstack_methods"));
		output_file = GTK_ENTRY(gtk_builder_get_object(gui.builder, "entryresultfile"));
		overwrite = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkbutoverwrite"));
		sigSpin[0] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
		sigSpin[1] = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));
		rejec_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comborejection"));
		norm_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combonormalize"));
		weighing_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboweighing"));
		force_norm = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkforcenorm"));
		fast_norm = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "checkfastnorm"));
		max_framing = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "check_maximize_framing"));
		norm_to_max = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_normalise_to_max"));
		RGB_equal = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_RGBequal"));
		rejmaps = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "rejmaps_checkbutton"));
		merge_rejmaps = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "merge_rejmaps_checkbutton"));
		upscale_at_stacking = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "check_upscale_at_stacking"));
		feather_dist = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_stack_feather_dist"));
		blend_frame = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_blend_frame"));
		overlap_norm = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "check_norm_overlap"));
		force32b = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "check_force32b"));
		/* STACK_MPP widgets */
		mpp_scale_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_mpp_scale"));
		mpp_engine_combo  = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combo_mpp_engine"));
		mpp_stack_percent = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_stack_percent"));
		mpp_stack_frames  = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_stack_frames"));
		mpp_bg_fraction   = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_bg_fraction"));
		mpp_bg_blend      = GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "spin_mpp_bg_blend"));
		mpp_skip_failed   = GTK_CHECK_BUTTON(gtk_builder_get_object(gui.builder, "check_mpp_skip_failed"));
	}

	if (processing_is_job_active()) {
		PRINT_ANOTHER_THREAD_RUNNING;
		return;
	}

	stackparam.sig[0] = (float) gtk_spin_button_get_value(sigSpin[0]);
	stackparam.sig[1] = (float) gtk_spin_button_get_value(sigSpin[1]);
	stackparam.type_of_rejection = gtk_drop_down_get_selected(rejec_combo);
	stackparam.create_rejmaps = siril_toggle_get_active(GTK_WIDGET(rejmaps));
	stackparam.merge_lowhigh_rejmaps = siril_toggle_get_active(GTK_WIDGET(merge_rejmaps));
	stackparam.normalize = gtk_drop_down_get_selected(norm_combo);
	stackparam.force_norm = siril_toggle_get_active(GTK_WIDGET(force_norm));
	stackparam.output_norm = siril_toggle_get_active(GTK_WIDGET(norm_to_max)) && gtk_widget_is_visible(norm_to_max);
	stackparam.maximize_framing = siril_toggle_get_active(GTK_WIDGET(max_framing)) && gtk_widget_is_visible(GTK_WIDGET(max_framing));
	stackparam.upscale_at_stacking = siril_toggle_get_active(GTK_WIDGET(upscale_at_stacking)) && gtk_widget_is_visible(GTK_WIDGET(upscale_at_stacking));
	stackparam.coeff.offset = NULL;
	stackparam.coeff.mul = NULL;
	stackparam.coeff.scale = NULL;
	stackparam.method =	stacking_methods[gtk_drop_down_get_selected(method_combo)];
	gboolean weighing_is_enabled = gtk_widget_get_visible(GTK_WIDGET(weighing_combo));
	if (weighing_is_enabled) {
		int weight_type = gtk_drop_down_get_selected(weighing_combo);
		if (weight_type == NOISE_WEIGHT) {
			int norm_type = gtk_drop_down_get_selected(norm_combo);
			if (norm_type == NO_NORM) {
				siril_log_error(_("Weighting by noise is allowed only if normalization has been activated, ignoring weights.\n"));
				weight_type = NO_WEIGHT;
			} else
				stackparam.weighting_type = norm_type;
		} else {
			stackparam.weighting_type = weight_type;
		}
	}
	stackparam.equalizeRGB = siril_toggle_get_active(GTK_WIDGET(RGB_equal)) && gtk_widget_is_visible(RGB_equal)  && (gtk_drop_down_get_selected(norm_combo) != NO_NORM);
	stackparam.lite_norm = siril_toggle_get_active(GTK_WIDGET(fast_norm));
	stackparam.feather_dist = (int)gtk_spin_button_get_value(feather_dist) * gtk_widget_get_visible(blend_frame);
	stackparam.overlap_norm = siril_toggle_get_active(GTK_WIDGET(overlap_norm)) * gtk_widget_get_visible(blend_frame);
	if (stackparam.overlap_norm && stackparam.weighting_type == NOISE_WEIGHT) {
		siril_log_error(_("Weighting by noise cannot be used with overlap normalization, ignoring weights.\n"));
		stackparam.weighting_type = NO_WEIGHT;
	}

	stackparam.use_32bit_output = siril_toggle_get_active(GTK_WIDGET(force32b)) || evaluate_stacking_should_output_32bits(stackparam.method,
			&com.seq, stackparam.nb_images_to_stack, &error);
	if (error) {
		siril_log_error(error);
		return;
	}

	// ensure we have no normalization if not supported by the stacking method
	if (stackparam.method != stack_median && stackparam.method != stack_mean_with_rejection)
		stackparam.normalize = NO_NORM;
	stackparam.seq = &com.seq;
	stackparam.reglayer = get_registration_layer(stackparam.seq);
	/* STACK_MPP doesn't go through the homography-aware regdata path — it
	 * reads its own per-AP per-frame shifts from the .mpp sidecar. The
	 * "registration data with more than simple shifts" prompt would fire
	 * for a sidecar's quality-only regdata so we skip it for STACK_MPP. */
	if (stackparam.method != stack_mpp_handler) {
		// checking regdata is absent, or if present, is only shift
		if (!test_regdata_is_valid_and_shift(stackparam.seq, stackparam.reglayer)) {
			int confirm = siril_confirm_dialog(_("Registration data found"),
				_("Stacking has detected registration data with more than simple shifts.\n"
				"Normally, you should apply existing registration before stacking."),
				_("Stack anyway"));
			if (!confirm)
				return;
		}
	}

	if (stackparam.overlap_norm && stackparam.nb_images_to_stack > 20) {
		gchar *onorm_msg = g_strdup_printf(_("You have chosen to compute normalization on overlaps with more than %d images.\n"
			"This option should normally be used to stitch stacked mosaic tiles, not subs.\n"
			"If you proceed, execution may be slow"), MAX_IMAGES_FOR_OVERLAP);
		int confirm = siril_confirm_dialog(_("Large number of images"),
			onorm_msg,
			_("Stack anyway"));
		g_free(onorm_msg);
		if (!confirm)
			return;
	}

	/* Do not display that cause it uses the generic function that already
	 * displays this text
	 */
	if (stackparam.method != &stack_summing_generic)
		siril_log_info(_("Stacking: processing...\n"));
	gettimeofday(&stackparam.t_start, NULL);
	set_cursor_waiting(TRUE);

	stackparam.output_overwrite = siril_toggle_get_active(GTK_WIDGET(overwrite));
	stackparam.output_filename = gtk_editable_get_text(GTK_EDITABLE(output_file));

	/* Stacking. Result is in gfit if success */
	struct stacking_args *params = calloc(1, sizeof(struct stacking_args));
	stacking_args_deep_copy(&stackparam, params);

	/* For STACK_MPP, capture the stack-side widget values into a fresh
	 * mpp_config_t and hang it off params. stack_mpp_handler reads from
	 * params->mpp_cfg and stacking_args_deep_free will release it.
	 * Allocated AFTER deep_copy so stackparam never owns the pointer. */
	if (stackparam.method == stack_mpp_handler) {
		mpp_config_t *cfg = calloc(1, sizeof(*cfg));
		mpp_config_defaults(cfg);
		const guint scale_idx = gtk_drop_down_get_selected(mpp_scale_combo);
		if (scale_idx < G_N_ELEMENTS(mpp_output_scales))
			cfg->output_scale = mpp_output_scales[scale_idx];
		else	/* GTK_INVALID_LIST_POSITION */
			cfg->output_scale = 1.0;
		cfg->stack_frame_percent                     = gtk_spin_button_get_value_as_int(mpp_stack_percent);
		cfg->stack_frame_number                      = gtk_spin_button_get_value_as_int(mpp_stack_frames);
		cfg->stack_frames_background_fraction        = gtk_spin_button_get_value(mpp_bg_fraction);
		cfg->stack_frames_background_blend_threshold = gtk_spin_button_get_value(mpp_bg_blend);
		cfg->stack_skip_failed_aps                   = siril_toggle_get_active(GTK_WIDGET(mpp_skip_failed));
		if (mpp_engine_combo) {
			/* Combo item indices match enum mpp_stack_method. */
			const guint eng = gtk_drop_down_get_selected(mpp_engine_combo);
			cfg->stack_method = (eng == 1) ? MPP_STACK_WARP : MPP_STACK_PATCH;
		}
		params->mpp_cfg = cfg;
	}

	if (!start_in_new_thread(stack_function_handler, params)) {
		stacking_args_deep_free(params);
	}
}

void on_seqstack_button_clicked (GtkButton *button, gpointer user_data){
	control_window_switch_to_tab(OUTPUT_LOGS);
	if (stack_lsf_append_mode())
		start_lsf_append();
	else if (stack_lsf_mode())
		start_lsf_stacking();
	else start_stacking();
}

void on_comboboxstack_methods_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	static GtkNotebook* notebook = NULL;
	if (!notebook)
		notebook = GTK_NOTEBOOK(gtk_builder_get_object(gui.builder, "notebook4"));
	com.pref.stack.method = gtk_drop_down_get_selected(box);

	gtk_notebook_set_current_page(notebook, com.pref.stack.method);
	update_stack_interface(TRUE);
}

void on_combo_mpp_engine_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)pspec;
	(void)user_data;
	static GtkWidget *skip_failed = NULL;
	if (!skip_failed)
		skip_failed = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_mpp_skip_failed"));
	if (!skip_failed)
		return;
	/* "Skip failed AP shifts" is a patch-engine policy: the warp engine
	 * handles failed pairs itself (smoothed shift predictions feed the
	 * displacement field) and never reads the flag. */
	const gboolean patch =
	    gtk_drop_down_get_selected(GTK_DROP_DOWN(obj)) != 1;
	gtk_widget_set_sensitive(skip_failed, patch);
}

void on_combonormalize_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	(void)obj;
	(void)pspec;
	static GtkWidget *widgetnormalize = NULL, *force_norm = NULL, *fast_norm = NULL;
	if (!widgetnormalize) {
		widgetnormalize = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combonormalize"));
		force_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkforcenorm"));
		fast_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkfastnorm"));
	}
	gtk_widget_set_sensitive(force_norm, gtk_drop_down_get_selected(GTK_DROP_DOWN(widgetnormalize)) != 0);
	gtk_widget_set_sensitive(fast_norm, gtk_drop_down_get_selected(GTK_DROP_DOWN(widgetnormalize)) != 0);
}

void on_comboweighing_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	com.pref.stack.weighting_method = gtk_drop_down_get_selected(box);
}

void on_stack_siglow_button_value_changed(GtkSpinButton *button, gpointer user_data) {
	rejection type_of_rejection = gtk_drop_down_get_selected((GtkDropDown *)user_data);
	if (stack_lsf_mode())
		update_lsf_hist_widgets();

	switch (type_of_rejection) {
	case PERCENTILE:
		com.pref.stack.percentile_low = gtk_spin_button_get_value(button);
		break;
	case LINEARFIT:
		com.pref.stack.linear_low = gtk_spin_button_get_value(button);
		break;
	case SIGMA:
	case MAD:
	case SIGMEDIAN:
	case WINSORIZED:
		com.pref.stack.sigma_low = gtk_spin_button_get_value(button);
		break;
	case GESDT:
		break;
	default:
		return;
	}
}

void on_stack_sighigh_button_value_changed(GtkSpinButton *button, gpointer user_data) {
	rejection type_of_rejection = gtk_drop_down_get_selected((GtkDropDown *)user_data);
	if (stack_lsf_mode())
		update_lsf_hist_widgets();

	switch (type_of_rejection) {
	case PERCENTILE:
		com.pref.stack.percentile_high = gtk_spin_button_get_value(button);
		break;
	case LINEARFIT:
		com.pref.stack.linear_high = gtk_spin_button_get_value(button);
		break;
	case SIGMA:
	case MAD:
	case SIGMEDIAN:
	case WINSORIZED:
		com.pref.stack.sigma_high = gtk_spin_button_get_value(button);
		break;
	case GESDT:
		break;
	default:
		return;
	}
}

void on_comborejection_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *box = GTK_DROP_DOWN(obj);
	(void)pspec;
	rejection type_of_rejection = gtk_drop_down_get_selected(box);
	static GtkWidget *labellow = NULL, *labelhigh = NULL;
	static GtkWidget *siglow = NULL, *sighigh = NULL;
	static GtkWidget *rejmaps = NULL, *merge_rejmaps = NULL;

	if (!labellow) {
		labellow = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_low"));
		labelhigh = GTK_WIDGET(gtk_builder_get_object(gui.builder, "label_high"));
		siglow = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
		sighigh = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));
		rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "rejmaps_checkbutton"));
		merge_rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "merge_rejmaps_checkbutton"));
	}

	g_signal_handlers_block_by_func(GTK_SPIN_BUTTON(siglow), on_stack_siglow_button_value_changed, NULL);
	g_signal_handlers_block_by_func(GTK_SPIN_BUTTON(sighigh), on_stack_sighigh_button_value_changed, NULL);
	/* set default values */
	switch (type_of_rejection) {
		case NO_REJEC:
			gtk_widget_set_visible(siglow, FALSE);
			gtk_widget_set_visible(rejmaps, FALSE);
			gtk_widget_set_visible(sighigh, FALSE);
			gtk_widget_set_visible(labellow, FALSE);
			gtk_widget_set_visible(labelhigh, FALSE);
			gtk_widget_set_visible(merge_rejmaps, FALSE);
			break;
		case PERCENTILE:
			gtk_widget_set_visible(siglow, TRUE);
			gtk_widget_set_visible(rejmaps, TRUE);
			gtk_widget_set_visible(sighigh, TRUE);
			gtk_widget_set_visible(labellow, TRUE);
			gtk_widget_set_visible(labelhigh, TRUE);
			gtk_widget_set_visible(merge_rejmaps, TRUE);
			gtk_widget_set_tooltip_text(siglow, _("Low clipping factor for the percentile clipping rejection algorithm."));
			gtk_widget_set_tooltip_text(sighigh, _("High clipping factor for the percentile clipping rejection algorithm."));
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(siglow), com.pref.stack.percentile_low);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(sighigh), com.pref.stack.percentile_high);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(siglow), 0.0, 1.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(sighigh), 0.0, 1.0);
			gtk_label_set_text(GTK_LABEL(labellow), _("Percentile low: "));
			gtk_label_set_text(GTK_LABEL(labelhigh), _("Percentile high: "));
			break;
		case LINEARFIT:
			gtk_widget_set_visible(siglow, TRUE);
			gtk_widget_set_visible(rejmaps, TRUE);
			gtk_widget_set_visible(sighigh, TRUE);
			gtk_widget_set_visible(labellow, TRUE);
			gtk_widget_set_visible(labelhigh, TRUE);
			gtk_widget_set_visible(merge_rejmaps, TRUE);
			gtk_widget_set_tooltip_text(siglow, _("Tolerance for low pixel values of the linear fit clipping algorithm, in sigma units."));
			gtk_widget_set_tooltip_text(sighigh, _("Tolerance for high pixel values of the linear fit clipping algorithm, in sigma units."));
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(siglow), 0.0, 10.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(sighigh), 0.0, 10.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(siglow), com.pref.stack.linear_low);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(sighigh), com.pref.stack.linear_high);
			gtk_label_set_text(GTK_LABEL(labellow), _("Linear low: "));
			gtk_label_set_text(GTK_LABEL(labelhigh), _("Linear high: "));
			break;
		default:
		case SIGMA:
		case MAD:
		case SIGMEDIAN:
		case WINSORIZED:
			gtk_widget_set_visible(siglow, TRUE);
			gtk_widget_set_visible(rejmaps, TRUE);
			gtk_widget_set_visible(sighigh, TRUE);
			gtk_widget_set_visible(labellow, TRUE);
			gtk_widget_set_visible(labelhigh, TRUE);
			gtk_widget_set_visible(merge_rejmaps, TRUE);
			gtk_widget_set_tooltip_text(siglow, _("Low clipping factor for the sigma clipping rejection algorithm."));
			gtk_widget_set_tooltip_text(sighigh, _("High clipping factor for the sigma clipping rejection algorithm."));
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(siglow), 0.0, 10.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(sighigh), 0.0, 10.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(siglow), com.pref.stack.sigma_low);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(sighigh), com.pref.stack.sigma_high);
			gtk_label_set_text(GTK_LABEL(labellow), _("Sigma low: "));
			gtk_label_set_text(GTK_LABEL(labelhigh), _("Sigma high: "));
			break;
		case GESDT:
			gtk_widget_set_visible(siglow, TRUE);
			gtk_widget_set_visible(rejmaps, TRUE);
			gtk_widget_set_visible(sighigh, TRUE);
			gtk_widget_set_visible(labellow, TRUE);
			gtk_widget_set_visible(labelhigh, TRUE);
			gtk_widget_set_visible(merge_rejmaps, TRUE);
			gtk_widget_set_tooltip_text(siglow, _("Expected fraction of maximum outliers for the Generalized Extreme Studentized Deviate Test algorithm."));
			gtk_widget_set_tooltip_text(sighigh, _("Probability of making a false positive for the Generalized Extreme Studentized Deviate Test algorithm. "
					"Increasing this value will reject more pixels."));
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(siglow), 0.0, 1.0);
			gtk_spin_button_set_range(GTK_SPIN_BUTTON(sighigh), 0.0, 1.0);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(siglow), 0.3);
			gtk_spin_button_set_value(GTK_SPIN_BUTTON(sighigh), 0.05);
			gtk_label_set_text(GTK_LABEL(labellow), _("ESD Outliers: "));
			gtk_label_set_text(GTK_LABEL(labelhigh), _("ESD Significance: "));
	}

	g_signal_handlers_unblock_by_func(GTK_SPIN_BUTTON(siglow), on_stack_siglow_button_value_changed, NULL);
	g_signal_handlers_unblock_by_func(GTK_SPIN_BUTTON(sighigh), on_stack_sighigh_button_value_changed, NULL);

	com.pref.stack.rej_method = gtk_drop_down_get_selected(box);
	/* a Living Stack has its own rejection maps option, stack_lsf_rejmaps */
	if (stack_lsf_mode()) {
		gtk_widget_set_visible(rejmaps, FALSE);
		gtk_widget_set_visible(merge_rejmaps, FALSE);
		update_lsf_hist_widgets();
	}
}

void on_rejmaps_toggled(GtkCheckButton *button, gpointer user_data) {
	static GtkWidget *merge = NULL;
	if (!merge) merge = GTK_WIDGET(gtk_builder_get_object(gui.builder, "merge_rejmaps_checkbutton"));
	gtk_widget_set_sensitive(merge, siril_toggle_get_active(GTK_WIDGET(button)));
}

void on_check_maximize_framing_toggled(GtkCheckButton *button, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_stacksel_changed(GObject *obj, GParamSpec *pspec, gpointer user_data) {
	GtkDropDown *widget = GTK_DROP_DOWN(obj);
	(void)pspec;
	const gchar *caller = gtk_buildable_get_buildable_id(GTK_BUILDABLE (widget));
	if (g_str_has_prefix(caller, "filter_type")) {
		update_adjustment = (int)g_ascii_strtod(caller + 11, NULL) - 1; // filter_type1, 2 or 3 to be parsed as 0, 1 or 2
	}
	update_stack_interface(TRUE);
}

void on_spinbut_percent_change(GtkSpinButton *spinbutton, gpointer user_data) {
	update_stack_interface(TRUE);
}

void on_filter_add1_clicked(GtkButton *button, gpointer user_data){
	static GtkWidget *combofilter2 = NULL, *stackspin2 = NULL, *filter_add2 = NULL;
	static GtkWidget *filter_rem2 = NULL, *labelfilter2 = NULL, *filter_type2 = NULL;
	if (!combofilter2) {
		combofilter2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combofilter2"));
		stackspin2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin2"));
		filter_add2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_add2"));
		filter_rem2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_rem2"));
		labelfilter2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelfilter2"));
		filter_type2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type2"));
	}
	gtk_widget_set_visible(combofilter2, TRUE);
	gtk_widget_set_visible(stackspin2, TRUE);
	gtk_widget_set_visible(filter_add2, TRUE);
	gtk_widget_set_visible(filter_rem2, TRUE);
	gtk_widget_set_visible(labelfilter2, TRUE);
	gtk_widget_set_visible(filter_type2, TRUE);
	update_stack_interface(TRUE);
}

void on_filter_add2_clicked(GtkButton *button, gpointer user_data){
	static GtkWidget *combofilter3 = NULL, *stackspin3 = NULL, *filter_rem3 = NULL;
	static GtkWidget *labelfilter3 = NULL, *filter_type3 = NULL;
	if (!combofilter3) {
		combofilter3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combofilter3"));
		stackspin3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin3"));
		filter_rem3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_rem3"));
		labelfilter3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelfilter3"));
		filter_type3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type3"));
	}
	gtk_widget_set_visible(combofilter3, TRUE);
	gtk_widget_set_visible(stackspin3, TRUE);
	gtk_widget_set_visible(filter_rem3, TRUE);
	gtk_widget_set_visible(labelfilter3, TRUE);
	gtk_widget_set_visible(filter_type3, TRUE);
	update_stack_interface(TRUE);
}

void on_filter_rem2_clicked(GtkButton *button, gpointer user_data){
	static GtkWidget *combofilter2 = NULL, *stackspin2 = NULL, *filter_add2 = NULL;
	static GtkWidget *filter_rem2 = NULL, *labelfilter2 = NULL, *filter_type2 = NULL;
	if (!combofilter2) {
		combofilter2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combofilter2"));
		stackspin2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin2"));
		filter_add2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_add2"));
		filter_rem2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_rem2"));
		labelfilter2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelfilter2"));
		filter_type2 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type2"));
	}
	gtk_widget_set_visible(combofilter2, FALSE);
	gtk_widget_set_visible(stackspin2, FALSE);
	gtk_widget_set_visible(filter_add2, FALSE);
	gtk_widget_set_visible(filter_rem2, FALSE);
	gtk_widget_set_visible(labelfilter2, FALSE);
	gtk_widget_set_visible(filter_type2, FALSE);
	update_stack_interface(TRUE);
}

void on_filter_rem3_clicked(GtkButton *button, gpointer user_data){
	static GtkWidget *combofilter3 = NULL, *stackspin3 = NULL, *filter_rem3 = NULL;
	static GtkWidget *labelfilter3 = NULL, *filter_type3 = NULL;
	if (!combofilter3) {
		combofilter3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combofilter3"));
		stackspin3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin3"));
		filter_rem3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_rem3"));
		labelfilter3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "labelfilter3"));
		filter_type3 = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type3"));
	}
	gtk_widget_set_visible(combofilter3, FALSE);
	gtk_widget_set_visible(stackspin3, FALSE);
	gtk_widget_set_visible(filter_rem3, FALSE);
	gtk_widget_set_visible(labelfilter3, FALSE);
	gtk_widget_set_visible(filter_type3, FALSE);
	update_stack_interface(TRUE);
}

void get_sequence_filtering_from_gui(seq_image_filter *filtering_criterion,
		double *filtering_parameter) {
	int filter, guifilter, channel = 0, type;
	gboolean is_ksig = FALSE;
	double percent = 0.0;
	static GtkDropDown *filter_combo[] = {NULL, NULL, NULL};
	static GtkAdjustment *stackadj[] = {NULL, NULL, NULL};
	static GtkWidget *spin[] = {NULL, NULL, NULL};
	static GtkWidget *ksig[] = {NULL, NULL, NULL};
	if (!spin[0]) {
		spin[0] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin1"));
		spin[1] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin2"));
		spin[2] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stackspin3"));
		stackadj[0] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[0]));
		stackadj[1] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[1]));
		stackadj[2] = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(spin[2]));
		filter_combo[0] = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combofilter1"));
		filter_combo[1] = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combofilter2"));
		filter_combo[2] = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combofilter3"));
		ksig[0] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type1"));
		ksig[1] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type2"));
		ksig[2] = GTK_WIDGET(gtk_builder_get_object(gui.builder, "filter_type3"));
	}
	for (filter = 0, guifilter = 0; guifilter < 3; guifilter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[guifilter]))) {
			continue;
		}

		type = gtk_drop_down_get_selected(filter_combo[guifilter]);
		if (type != ALL_IMAGES && type != SELECTED_IMAGES) {
			channel = get_registration_layer(&com.seq);
			percent = gtk_adjustment_get_value(stackadj[guifilter]);
			is_ksig = gtk_drop_down_get_selected(GTK_DROP_DOWN(ksig[guifilter]));
		}
		if (update_adjustment == filter) {
			g_signal_handlers_block_by_func(stackadj[filter], on_stacksel_changed, NULL);
			gtk_adjustment_set_upper(stackadj[filter], filter_maxvals[is_ksig]);
			gtk_adjustment_set_value(stackadj[filter], filter_initvals[is_ksig]);
			gtk_adjustment_set_step_increment(stackadj[filter], filter_increments[is_ksig]);
			gtk_widget_set_tooltip_text(spin[filter], filter_tooltip_text[is_ksig]);
			g_signal_handlers_unblock_by_func(stackadj[filter], on_stacksel_changed, NULL);
			update_adjustment = -1;
		}
		switch (type) {
			default:
			case ALL_IMAGES:
				stackfilters[filter].filter = seq_filter_all;
				stackfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				gtk_widget_set_visible(ksig[guifilter], FALSE);
				break;
			case SELECTED_IMAGES:
				stackfilters[filter].filter = seq_filter_included;
				stackfilters[filter].param = 0.0;
				gtk_widget_set_visible(spin[guifilter], FALSE);
				gtk_widget_set_visible(ksig[guifilter], FALSE);
				break;
			case BEST_PSF_IMAGES:
				stackfilters[filter].filter = seq_filter_fwhm;
				stackfilters[filter].param = compute_highest_accepted_fwhm(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_WPSF_IMAGES:
				stackfilters[filter].filter = seq_filter_weighted_fwhm;
				stackfilters[filter].param = compute_highest_accepted_weighted_fwhm(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_ROUND_IMAGES:
				stackfilters[filter].filter = seq_filter_roundness;
				stackfilters[filter].param = compute_lowest_accepted_roundness(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_BKG_IMAGES:
				stackfilters[filter].filter = seq_filter_background;
				stackfilters[filter].param = compute_highest_accepted_background(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_NBSTARS_IMAGES:
				stackfilters[filter].filter = seq_filter_nbstars;
				stackfilters[filter].param = compute_lowest_accepted_nbstars(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
			case BEST_QUALITY_IMAGES:
				stackfilters[filter].filter = seq_filter_quality;
				stackfilters[filter].param = compute_lowest_accepted_quality(
						stackparam.seq, channel, percent, is_ksig);
				gtk_widget_set_visible(spin[guifilter], TRUE);
				gtk_widget_set_visible(ksig[guifilter], TRUE);
				break;
		}
		filter++;
	}
	stackfilters[filter].filter = NULL;

	if (filter == 1) {
		*filtering_criterion = stackfilters[0].filter;
		*filtering_parameter = stackfilters[0].param;
	} else {
		*filtering_criterion = create_multiple_filter_from_list(stackfilters);
		*filtering_parameter = 0.0;
	}
}

static void update_filter_label() {
	static GtkDropDown *filter_combo[3] = { NULL };
	static GtkLabel *filter_label[3] = { NULL };
	if (!filter_combo[0]) {
		filter_combo[0] = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combofilter1"));
		filter_combo[1] = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combofilter2"));
		filter_combo[2] = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combofilter3"));
		filter_label[0] = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelfilter1"));
		filter_label[1] = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelfilter2"));
		filter_label[2] = GTK_LABEL(gtk_builder_get_object(gui.builder, "labelfilter3"));
	}

	for (int filter = 0; filter < 3; filter++) {
		if (!gtk_widget_get_visible(GTK_WIDGET(filter_combo[filter]))) {
			break;
		}

		int type = gtk_drop_down_get_selected(filter_combo[filter]);
		double param = stackfilters[filter].param;
		gchar *filter_str;
		if (param == DBL_MIN || param == DBL_MAX || param == 0.0) {
			if (type == ALL_IMAGES || type == SELECTED_IMAGES)
				filter_str = g_strdup("");
			else filter_str = g_strdup("N/A");
		} else {
			switch (type) {
			default:
			case ALL_IMAGES:
			case SELECTED_IMAGES:
				filter_str = g_strdup("");
				break;
			case BEST_PSF_IMAGES:
			case BEST_WPSF_IMAGES:
				filter_str = g_strdup_printf("< %.2lf", param);
				break;
			case BEST_BKG_IMAGES :
				filter_str = (param < 1.) ? g_strdup_printf("< %.5lf", param) : g_strdup_printf("< %d", (int)param);
				break;
			case BEST_ROUND_IMAGES:
			case BEST_QUALITY_IMAGES:
				filter_str = g_strdup_printf("> %.3lf", param);
				break;
			case BEST_NBSTARS_IMAGES:
				filter_str = g_strdup_printf("> %d", (int)param);
				break;
			}
		}
		gtk_label_set_text(filter_label[filter], filter_str);
		g_free(filter_str);
	}
}

/* Activates or not the stack button if there are 2 or more selected images,
 * all data related to stacking is set in stackparam, except the method itself,
 * determined at stacking start.
 */
void update_stack_interface(gboolean dont_change_stack_type) {
	static GtkWidget *go_stack = NULL, *widgetnormalize = NULL, *force_norm =
			NULL, *output_norm = NULL, *RGB_equal = NULL, *fast_norm = NULL, *max_framing = NULL,
			*upscale_at_stacking = NULL, *blend_frame = NULL, *overlap_norm = NULL;
	static GtkWidget *lsf_grid = NULL, *lsf_framing = NULL, *force32b = NULL, *rejmaps = NULL,
			*merge_rejmaps = NULL;
	static GtkWidget *append_box = NULL, *output_box = NULL, *lsf_framing_label = NULL, *lsf_comp_label = NULL,
			*lsf_comp = NULL, *lsf_hist_same = NULL, *lsf_rejmaps = NULL, *weighing_combo = NULL,
			*siglow = NULL, *sighigh = NULL, *append_framing = NULL, *append_noreg = NULL,
			*expand_step = NULL, *expand_force = NULL, *info_label = NULL, *force_reg = NULL;
	static GtkDropDown *rejec_combo = NULL;
	static GtkDropDown *method_combo = NULL, *filter_combo = NULL;
	static GtkLabel *result_label = NULL;
	static GtkExpander *stack_expander_method = NULL, *stack_expander_output = NULL;
	static GtkWidget *seq_filters_frame = NULL;
	gchar *labelbuffer;

	if(!go_stack) {
		go_stack = GTK_WIDGET(gtk_builder_get_object(gui.builder, "gostack_button"));
		filter_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combofilter1"));
		seq_filters_frame = GTK_WIDGET(gtk_builder_get_object(gui.builder, "seq_filters_frame_stack"));
		method_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboboxstack_methods"));
		widgetnormalize = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combonormalize"));
		force_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkforcenorm"));
		fast_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkfastnorm"));
		result_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "stackfilter_label"));
		output_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_normalise_to_max"));
		RGB_equal = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_RGBequal"));
		max_framing = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_maximize_framing"));
		upscale_at_stacking = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_upscale_at_stacking"));
		blend_frame = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_blend_frame"));
		overlap_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_norm_overlap"));
		stack_expander_method = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "stack_expander_method"));
		stack_expander_output = GTK_EXPANDER(gtk_builder_get_object(gui.builder, "stack_expander_output"));
		lsf_grid = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_grid"));
		lsf_framing = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_framing"));
		force32b = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_force32b"));
		rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "rejmaps_checkbutton"));
		merge_rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "merge_rejmaps_checkbutton"));
		rejec_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comborejection"));
		append_box = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_box"));
		output_box = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_output_box"));
		lsf_framing_label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_framing_label"));
		lsf_comp_label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_compression_label"));
		lsf_comp = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_compression"));
		lsf_hist_same = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_hist_same"));
		lsf_rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps"));
		weighing_combo = GTK_WIDGET(gtk_builder_get_object(gui.builder, "comboweighing"));
		siglow = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_siglow_button"));
		sighigh = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_sighigh_button"));
		append_framing = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_framing"));
		append_noreg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_noreg"));
		expand_step = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_expand_step"));
		expand_force = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_expand_force"));
		info_label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_info_label"));
		force_reg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_force_reg"));
	}
	gboolean seqloaded = sequence_is_loaded();
	gtk_widget_set_sensitive(GTK_WIDGET(stack_expander_method), seqloaded);
	gtk_expander_set_expanded(stack_expander_method, seqloaded);
	gtk_widget_set_sensitive(GTK_WIDGET(stack_expander_output), seqloaded);
	gtk_expander_set_expanded(stack_expander_output, seqloaded);
	if (seqloaded) {
		expander_clear_stuck_insensitive(GTK_WIDGET(stack_expander_method));
		expander_clear_stuck_insensitive(GTK_WIDGET(stack_expander_output));
	}
	if (!seqloaded) {
		return;
	}
	stackparam.seq = &com.seq;

	if (!dont_change_stack_type && stackparam.seq->selnum < stackparam.seq->number) {
		g_signal_handlers_block_by_func(filter_combo, on_stacksel_changed, NULL);
		gtk_drop_down_set_selected(filter_combo, SELECTED_IMAGES);
		g_signal_handlers_unblock_by_func(filter_combo, on_stacksel_changed, NULL);
	}

	/* If a complete (registered) .mpp sidecar is loaded for this sequence
	 * (left by a previous "Multipoint Registration" run, or by the CLI
	 * register_mpp command), default the stack-method combo to STACK_MPP —
	 * the only method that can consume the sidecar. An analysis-only sidecar
	 * (Analyze / Stage A) carries no per-AP shifts and cannot be stacked, so
	 * don't steer the user into it; the cached run reflects whichever the
	 * sidecar held when the sequence was opened. dont_change_stack_type
	 * guards the case where the user has already manually picked something. */
	if (!dont_change_stack_type) {
		const mpp_run_t *run = mpp_get_cached_run();
		if (run && run->shifts)
			gtk_drop_down_set_selected(method_combo, STACK_MPP);
	}
	gboolean can_reframe = layer_has_usable_registration(&com.seq, get_registration_layer(&com.seq));
	gboolean can_upscale = can_reframe && !com.seq.is_variable && !com.seq.is_drizzle;
	gboolean must_reframe = can_reframe && com.seq.is_variable;

	int stack_method = gtk_drop_down_get_selected(method_combo);
	switch (stack_method) {
	default:
	case STACK_SUM:
	case STACK_MAX:
	case STACK_MIN:
		gtk_widget_set_sensitive(widgetnormalize, FALSE);
		gtk_widget_set_sensitive(force_norm, FALSE);
		gtk_widget_set_sensitive(fast_norm, FALSE);
		gtk_widget_set_visible(output_norm, FALSE);
		gtk_widget_set_visible(RGB_equal, FALSE);
		gtk_widget_set_visible(max_framing, can_reframe); // only shown if applicable
		if (can_reframe) {
			gtk_widget_set_sensitive(max_framing, !must_reframe);
			if (must_reframe)
				siril_toggle_set_active(GTK_WIDGET(max_framing), TRUE);
		}
		gtk_widget_set_visible(upscale_at_stacking, can_upscale); // only shown if applicable
		break;
	case STACK_MEAN:
	case STACK_MEDIAN:
		gtk_widget_set_sensitive(widgetnormalize, TRUE);
		gtk_widget_set_sensitive(force_norm,
				gtk_drop_down_get_selected(GTK_DROP_DOWN(widgetnormalize)) != 0);
		gtk_widget_set_sensitive(fast_norm,
				gtk_drop_down_get_selected(GTK_DROP_DOWN(widgetnormalize)) != 0);
		gtk_widget_set_visible(output_norm, TRUE);
		gtk_widget_set_visible(RGB_equal, TRUE);
		gtk_widget_set_visible(max_framing, can_reframe && stack_method != STACK_MEDIAN); // only shown if applicable and not for median
		gtk_widget_set_visible(upscale_at_stacking, can_upscale && stack_method != STACK_MEDIAN); // only shown if applicable and not for median
		if (can_reframe && stack_method != STACK_MEDIAN) {
			gtk_widget_set_sensitive(max_framing, !must_reframe);
			if (must_reframe)
				siril_toggle_set_active(GTK_WIDGET(max_framing), TRUE);
		}
	}

	gboolean can_feather = can_reframe && stack_method == STACK_MEAN;
	gtk_widget_set_visible(blend_frame, can_feather);
	if (can_feather) {
		gtk_widget_set_sensitive(overlap_norm, siril_toggle_get_active(GTK_WIDGET(max_framing)));
		if (!gtk_widget_get_sensitive(overlap_norm))
			siril_toggle_set_active(GTK_WIDGET(overlap_norm), FALSE);
	}

	/* A Living Stack is a mean with rejection made by lsf_init, which applies
	 * the registration itself with its own framing, always writes float32,
	 * and has its own rejection maps option (stack_lsf_rejmaps), as its maps
	 * are counts stored in the file, not images */
	gboolean lsf = stack_lsf_mode(), append = stack_lsf_append_mode();
	gtk_widget_set_visible(lsf_grid, lsf);
	gtk_widget_set_sensitive(GTK_WIDGET(method_combo), !lsf);
	gtk_button_set_label(GTK_BUTTON(go_stack), append ? _("Append to Living Stack") :
			lsf ? _("Create Living Stack") : _("Start stacking"));
	gtk_widget_set_visible(force32b, !lsf);
	gtk_widget_set_visible(force_norm, !lsf);
	/* an append takes these from the file (fill_widgets_from_lsf) */
	gtk_widget_set_visible(append_box, append);
	gtk_widget_set_visible(output_box, !append);
	gtk_widget_set_visible(lsf_framing_label, !append);
	gtk_widget_set_visible(lsf_framing, !append);
	gtk_widget_set_visible(lsf_comp_label, !append);
	gtk_widget_set_visible(lsf_comp, !append);
	gtk_widget_set_visible(lsf_hist_same, !append);
	/* facts of the file are locked; the rejection may be overridden (UI plan, C2) */
	GtkWidget *locked[] = { RGB_equal, output_norm, weighing_combo };
	for (int i = 0; i < (int) G_N_ELEMENTS(locked); i++)
		gtk_widget_set_sensitive(locked[i], !append);
	gtk_widget_set_sensitive(GTK_WIDGET(rejec_combo), TRUE);
	gtk_widget_set_sensitive(siglow, TRUE);
	gtk_widget_set_sensitive(sighigh, TRUE);
	if (append) {
		gtk_widget_set_sensitive(widgetnormalize, FALSE);
		gtk_widget_set_sensitive(fast_norm, FALSE);
	}
	if (lsf) {
		if (stack_method != STACK_MEAN) {
			/* the method handler comes back here */
			gtk_drop_down_set_selected(method_combo, STACK_MEAN);
			return;
		}
		gtk_widget_set_visible(max_framing, FALSE);
		gtk_widget_set_visible(upscale_at_stacking, FALSE);
		gtk_widget_set_visible(blend_frame, FALSE);
		gtk_widget_set_visible(rejmaps, FALSE);
		gtk_widget_set_visible(merge_rejmaps, FALSE);
		/* with the registration already applied, the frames define the grid */
		gtk_widget_set_sensitive(lsf_framing, can_reframe);
		update_lsf_hist_widgets();
		/* an append can only add rejection maps to a file that has none */
		gtk_widget_set_visible(lsf_rejmaps, !append || (append_info && !append_info->file.settings.rejmaps));
		gtk_check_button_set_label(GTK_CHECK_BUTTON(lsf_rejmaps), append ? _("Add rejection maps to the file, counted from this append") :
				_("Keep rejection maps in the file"));
	} else {
		/* as on_comborejection_changed shows them */
		gboolean has_rej = gtk_drop_down_get_selected(rejec_combo) != NO_REJEC;
		gtk_widget_set_visible(rejmaps, has_rej);
		gtk_widget_set_visible(merge_rejmaps, has_rej);
	}
	/* Multipoint stacking always honours the sequence's frame-selector
	 * inclusion state directly (in mpp_stack_apply); the global Stack-tab
	 * frame filter (and its frame-count label) is completely ignored here,
	 * so hide the whole section for STACK_MPP rather than leaving a stale
	 * "all" combo visible. A Living Stack uses the filters as stacking does
	 * (UI plan, A4): their list of frames is passed to the LSF job. */
	gtk_widget_set_visible(seq_filters_frame, stack_method != STACK_MPP);

	if (com.seq.reference_image == -1)
		com.seq.reference_image = sequence_find_refimage(&com.seq);
	stackparam.ref_image = com.seq.reference_image;

	get_sequence_filtering_from_gui(
			&stackparam.filtering_criterion, &stackparam.filtering_parameter);

	if (stackparam.description)
		g_free(stackparam.description);
	stackparam.description = describe_filter(stackparam.seq,
			stackparam.filtering_criterion, stackparam.filtering_parameter);

	update_filter_label();

	stackparam.nb_images_to_stack = compute_nb_filtered_images(&com.seq,
			stackparam.filtering_criterion, stackparam.filtering_parameter);
	/* an append may add a single frame */
	gboolean listed = stackparam.nb_images_to_stack >= (append ? 1 : 2) &&
		!stack_fill_list_of_unfiltered_images(&stackparam);

	if (append) {
		gboolean grow = gtk_drop_down_get_selected(GTK_DROP_DOWN(append_framing)) == 0;
		gboolean noreg = siril_toggle_get_active(append_noreg);
		gtk_widget_set_sensitive(append_framing, !noreg);
		gtk_widget_set_sensitive(expand_step, grow && !noreg);
		gtk_widget_set_sensitive(expand_force, grow && !noreg);
		gboolean fatal = TRUE, reusable = FALSE;
		gchar *problems = append_problems(&fatal);
		gchar *note = noreg ? NULL : append_registration_note(&reusable);
		gtk_widget_set_visible(force_reg, reusable);
		GString *text = g_string_new(append_summary ? append_summary : _("Choose a Living Stack file."));
		if (problems)
			g_string_append_printf(text, "\n%s", problems);
		if (note)
			g_string_append_printf(text, "\n%s", note);
		gtk_label_set_text(GTK_LABEL(info_label), text->str);
		g_string_free(text, TRUE);
		g_free(problems);
		g_free(note);
		labelbuffer = g_strdup_printf(_("Appending %d images of the %d of the sequence"),
				stackparam.nb_images_to_stack, com.seq.number);
		gtk_label_set_text(result_label, labelbuffer);
		g_free(labelbuffer);
		gtk_widget_set_sensitive(go_stack, append_info && !fatal && listed);
		return;
	}

	labelbuffer = g_strdup_printf(lsf ? _("Living Stack of %d images of the %d of the sequence") :
			_("Stacking %d images of the %d of the sequence"), stackparam.nb_images_to_stack, com.seq.number);
	gtk_label_set_text(result_label, labelbuffer);
	g_free(labelbuffer);
	gtk_widget_set_sensitive(go_stack, listed);
}
