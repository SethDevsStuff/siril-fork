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

/* The Living Stack modes of the Stacking tab: creating a Living Stack with
 * lsf_init, and appending to one with lsf_append. stacking.c calls the
 * stack_lsf_* and lsf_stacking_* functions. */

#include <math.h>
#include <string.h>

#include "algos/siril_wcs.h"
#include "core/gui_iface.h"
#include "core/processing.h"
#include "core/proto.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/gui_state.h"
#include "gui-gtk4/lsf_manager.h"
#include "gui-gtk4/lsf_stacking.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/stacking.h"
#include "gui-gtk4/utils.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "lsf/lsf.h"
#include "lsf/lsf_commands.h"
#include "lsf/lsf_fits.h"
#include "lsf/lsf_settings.h"
#include "lsf/lsf_weights.h"
#include "registration/registration.h"
#include "stacking/stacking.h"

/* Stacking result: an image, or a Living Stack (LSF) file made by lsf_init
 * from the same widgets, with the LSF options of stack_lsf_grid, or an
 * append to one with lsf_append, which takes its settings from the file */
enum { STACK_RESULT_IMAGE, STACK_RESULT_LSF, STACK_RESULT_LSF_IMAGE, STACK_RESULT_APPEND };

/* stack_lsf_framing items, in order */
static const framing_type lsf_framings[] = { FRAMING_MAX, FRAMING_CURRENT, FRAMING_MIN, FRAMING_COG };

/* The "Save in" entry holds a name pattern in both modes, parsed when the
 * result is written (path_parse). The LSF one starts from the stack-name
 * preference, as a stacking result does, with .lsf.fits. */
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
gboolean stack_lsf_mode(void) {
	return stack_result_mode() != STACK_RESULT_IMAGE;
}

gboolean stack_lsf_append_mode(void) {
	return stack_result_mode() == STACK_RESULT_APPEND;
}

/* the file chosen for an append, read by lsf_read_info, or NULL */
static lsf_info *append_info = NULL;
static guint append_serial = 0;		// changes each time append_info is read again
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
	/* a new file keeps rejection maps by default; an append
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

void lsf_stacking_update_hist_widgets(void) {
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
	lsf_stacking_update_hist_widgets();
}

/* the weight floor and range multiplier of a new file (Section 10.3): for
 * the weightings that freeze a range, wfwhm and nbstars. An append uses the
 * file's, which the summary shows. */
void lsf_stacking_update_weight_widgets(void) {
	static const char *ids[] = { "stack_lsf_weight_label", "stack_lsf_weight_floor", "stack_lsf_weight_range_label",
		"stack_lsf_weight_range" };
	weightingType w = gtk_drop_down_get_selected(GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboweighing")));
	gboolean shown = stack_lsf_mode() && !stack_lsf_append_mode() && (w == WFWHM_WEIGHT || w == NBSTARS_WEIGHT);
	for (int i = 0; i < (int) G_N_ELEMENTS(ids); i++)
		gtk_widget_set_visible(GTK_WIDGET(gtk_builder_get_object(gui.builder, ids[i])), shown);
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

/* blocks or unblocks the handlers of signal ("value-changed",
 * "notify::selected") on obj */
static void block_handlers(GObject *obj, const char *signal, gboolean block) {
	guint id;
	GQuark detail;
	if (!g_signal_parse_name(signal, G_OBJECT_TYPE(obj), &id, &detail, TRUE))
		return;
	GSignalMatchType match = G_SIGNAL_MATCH_ID | (detail ? G_SIGNAL_MATCH_DETAIL : 0);
	if (block)
		g_signal_handlers_block_matched(obj, match, id, detail, NULL, NULL, NULL);
	else g_signal_handlers_unblock_matched(obj, match, id, detail, NULL, NULL, NULL);
}

/* the widgets of the method expander show the settings of the file: the
 * facts (normalization, weighting, output normalization) cannot change, and
 * the rejection may be overridden. The handlers of the sigma values and of
 * the weighting only store the user's stacking preferences, so they are
 * blocked; the rejection handler sets the range of the sigma values and
 * runs, its preference is put back. */
static void fill_widgets_from_lsf(const lsf_settings *s) {
	GObject *siglow = gtk_builder_get_object(gui.builder, "stack_siglow_button");
	GObject *sighigh = gtk_builder_get_object(gui.builder, "stack_sighigh_button");
	GObject *weighing = gtk_builder_get_object(gui.builder, "comboweighing");
	int user_rej = com.pref.stack.rej_method;
	gtk_drop_down_set_selected(GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comborejection")), s->rej);
	com.pref.stack.rej_method = user_rej;
	block_handlers(siglow, "value-changed", TRUE);
	block_handlers(sighigh, "value-changed", TRUE);
	block_handlers(weighing, "notify::selected", TRUE);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(siglow), s->sig[0]);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(sighigh), s->sig[1]);
	gtk_drop_down_set_selected(GTK_DROP_DOWN(weighing), s->weighting);
	block_handlers(siglow, "value-changed", FALSE);
	block_handlers(sighigh, "value-changed", FALSE);
	block_handlers(weighing, "notify::selected", FALSE);
	gtk_drop_down_set_selected(GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "combonormalize")), s->norm);
	siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkfastnorm")), s->lite_norm);
	siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_RGBequal")), s->rgb_equal);
	siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_normalise_to_max")), s->outnorm);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_low")), s->hist_low);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_high")), s->hist_high);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_hist_margin")), s->hist_margin);
	siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps")), FALSE);
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
	append_serial++;
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
				/* a switch is chosen for one file */
				siril_toggle_set_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_weight_inverse")),
						FALSE);
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

/* after a registration of the sequence against a Living Stack:
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

/* stack_lsf_append_regmethod: how lsf_append registers (Section 10.4) */
enum { APPEND_REG_STARS, APPEND_REG_WCS, APPEND_REG_NONE };

static guint append_regmethod() {
	return gtk_drop_down_get_selected(GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_lsf_append_regmethod")));
}

/* the feathering distance of the Living Stack options, in the units of
 * stacking's Borders feathering, and their normalization on the overlaps */
static int blend_feather() {
	return gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder, "stack_lsf_feather")));
}

static gboolean blend_overlap_norm() {
	return siril_toggle_get_active(GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_overlap_norm")));
}

/* the masks of Siril's feathering are at a tenth of the image size
 * (MASK_SCALE): the distance in pixels the value gives */
static void on_lsf_feather_changed(GtkSpinButton *spin, gpointer label) {
	gchar *txt = g_strdup_printf(_("≈ %d px"), 10 * gtk_spin_button_get_value_as_int(spin));
	gtk_label_set_text(GTK_LABEL(label), txt);
	g_free(txt);
}

/* registration through the WCS needs the Stack's WCS, and plate-solved
 * frames, which the reference frame shows. Returns a message or NULL. */
static gchar *wcsreg_problems(gboolean *fatal) {
	*fatal = FALSE;
	if (!has_wcs(&append_info->file.primary)) {
		*fatal = TRUE;
		return g_strdup(_("This Living Stack has no WCS, so frames cannot be registered through theirs: open the file, "
					"plate-solve its stack (Tools > Astrometry) and write the solution into the file, or register "
					"with stars."));
	}
	fits meta = { 0 };
	gchar *msg = NULL;
	if (!seq_read_frame_metadata(&com.seq, sequence_find_refimage(&com.seq), &meta) && !has_wcs(&meta))
		msg = g_strdup(_("The reference frame of the sequence is not plate-solved: frames registered through their WCS "
					"must be. Frames that are not are not merged."));
	clearfits(&meta);
	return msg;
}

/* the frames to append beyond the weighting range frozen in the file
 * (Section 10.3), from the registration data of the sequence: an estimate,
 * since the append measures them again on the grid, against the file's
 * reference stars, and the log gives its own count. Returns a message, or
 * NULL when the file's weighting has no range or every frame is within it. */
static gchar *weight_range_note(void) {
	const lsf_settings *s = &append_info->file.settings;
	if ((s->weighting != WFWHM_WEIGHT && s->weighting != NBSTARS_WEIGHT) || s->weight_model != LSF_WMODEL_RANGE)
		return NULL;
	const struct stacking_args *sa = stack_gui_args();
	int layer = get_registration_layer(&com.seq);
	if (layer < 0 || !com.seq.regparam || !com.seq.regparam[layer] || !sa->image_indices || sa->nb_images_to_stack < 1)
		return g_strdup(_("The append measures the frames against the weighting range of the file; the log says "
					"if some are beyond it."));
	int n = sa->nb_images_to_stack;
	float *wfwhm = malloc(n * sizeof(float));
	int *nstars = malloc(n * sizeof(int));
	if (!wfwhm || !nstars) {
		free(wfwhm);
		free(nstars);
		return NULL;
	}
	for (int k = 0; k < n; k++) {
		const regdata *r = &com.seq.regparam[layer][sa->image_indices[k]];
		wfwhm[k] = r->weighted_fwhm;
		nstars[k] = r->number_of_stars;
	}
	int better = 0, worse = 0, floored = 0;
	lsf_weight_range_check(s, n, wfwhm, nstars, &better, &worse, &floored);
	free(wfwhm);
	free(nstars);
	if (better + worse == 0)
		return NULL;
	return g_strdup_printf(_("From the registration data of the sequence, about %d of the %d frames are beyond the "
				"weighting range of this Living Stack: %d better than the best of its first batch, which would "
				"all get the top weight, and %d worse than its worst, %d of them at the weight floor. Their "
				"weights would not follow their quality: switching to inverse weighting avoids it. The append "
				"measures the frames again, so its own count may differ."), better + worse, n, better, worse,
			floored);
}

/* the problems update_stack_interface shows and start_lsf_append checks:
 * the file, the batch and the pre-flight. Returns a message or NULL. */
static gchar *append_problems(gboolean *fatal) {
	*fatal = TRUE;
	if (!append_info)
		return NULL;
	if (append_info->no_append)
		return g_strdup_printf(_("This version cannot append to this file: %s."), append_info->no_append);
	gboolean batch_fatal = FALSE, flight_fatal = FALSE, wcs_fatal = FALSE;
	gchar *batch = lsf_check_batch(append_info, &com.seq, &batch_fatal);
	if (batch_fatal)
		return batch;
	/* a file made from an astrometric registration: star alignment cannot
	 * undistort with the frames' own SIP terms, only with the master */
	const lsf_settings *s = &append_info->file.settings;
	GtkDropDown *disto_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_lsf_append_disto"));
	if (s->undistorted && s->disto.index == DISTO_FILES && append_regmethod() == APPEND_REG_STARS
			&& gtk_drop_down_get_selected(disto_combo) != 1) {
		g_free(batch);
		return g_strdup(_("This Living Stack undistorts each frame with the SIP terms of its own WCS: plate-solve "
					"the frames and register them through the WCS, or undistort them with the distortion master."));
	}
	gchar *wcs = append_regmethod() == APPEND_REG_WCS ? wcsreg_problems(&wcs_fatal) : NULL;
	if (wcs_fatal) {
		g_free(batch);
		return wcs;
	}
	if (wcs) {
		gchar *both = batch ? g_strdup_printf("%s %s", batch, wcs) : g_strdup(wcs);
		g_free(batch);
		g_free(wcs);
		batch = both;
	}
	gchar *range = weight_range_note();
	if (range) {
		gchar *both = batch ? g_strdup_printf("%s %s", batch, range) : g_strdup(range);
		g_free(batch);
		g_free(range);
		batch = both;
	}
	GtkWidget *rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps"));
	gchar *out = append_out_path();
	gchar *flight = lsf_append_preflight(append_info, com.wd, out, stack_gui_args()->nb_images_to_stack,
			append_regmethod() == APPEND_REG_NONE, gtk_widget_get_visible(rejmaps) && siril_toggle_get_active(rejmaps),
			&flight_fatal);
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
 * tab or with register -extref=: what the append does with it.
 * Returns a note for the info label or NULL; *reusable when it is used. */
static gchar *append_registration_note(gboolean *reusable) {
	*reusable = FALSE;
	if (!append_info)
		return NULL;
	int layer;
	gchar *date = NULL, *why = NULL, *msg = NULL;
	const lsf_settings *s = &append_info->file.settings;
	*reusable = lsf_binding_check(&com.seq, &append_info->file, s->undistorted ? &s->disto : NULL, &layer, &date, &why);
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

/* the frames the stacking filters chose, for the LSF job, or
 * FALSE when there are none */
static gboolean set_lsf_frames(struct lsf_cmd_args *a) {
	if (!stack_gui_args()->image_indices || stack_gui_args()->nb_images_to_stack < 1)
		return FALSE;
	a->nb_images = stack_gui_args()->nb_images_to_stack;
	a->image_indices = g_memdup2(stack_gui_args()->image_indices, a->nb_images * sizeof(int));
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
	GtkWidget *rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps"));
	GtkDropDown *disto_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_lsf_append_disto"));
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

	/* overrides of the rejection and of the history test
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
	a->noreg = append_regmethod() == APPEND_REG_NONE;
	a->wcsreg = append_regmethod() == APPEND_REG_WCS;
	/* shown when a registration against the file can be reused */
	GtkWidget *force_reg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_force_reg"));
	GtkWidget *inverse = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_weight_inverse"));
	a->weight_inverse = gtk_widget_get_visible(inverse) && siril_toggle_get_active(inverse);
	a->force_reg = append_regmethod() == APPEND_REG_STARS && gtk_widget_get_visible(force_reg)
		&& siril_toggle_get_active(force_reg);
	/* an undistorted file: its model, or the master of the preferences */
	if (append_info->file.settings.undistorted && !a->noreg && gtk_drop_down_get_selected(disto_combo) == 1) {
		if (!com.pref.prepro.disto_lib || !com.pref.prepro.disto_lib[0]) {
			siril_message_dialog(GTK_MESSAGE_ERROR, _("No distortion master"), _("Set a distortion master in the "
						"preferences (Preprocessing), or use the model of the file."));
			lsf_cmd_args_free(a);
			return;
		}
		a->disto.index = DISTO_MASTER;
		a->disto.filename = g_strdup(com.pref.prepro.disto_lib);
	}
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
	/* feathering, and normalization on the overlap with the Stack */
	a->feather = blend_feather();
	a->overlap_norm = blend_overlap_norm();
	/* stack_lsf_append_compression: keep, none, gzip2, gzip1 */
	guint comp = gtk_drop_down_get_selected(comp_combo);
	if (comp >= 1 && comp <= 3) {
		a->compression = comp == 1 ? LSF_COMP_NONE : comp == 2 ? LSF_COMP_GZIP2 : LSF_COMP_GZIP1;
		a->compression_set = TRUE;
	}
	/* the plain image of the updated stack, named like the file written */
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

	/* memory and disk first */
	gboolean fatal = FALSE;
	gchar *flight = lsf_init_preflight(&com.seq, com.wd, path, stack_gui_args()->nb_images_to_stack,
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

	if (stack_gui_args()->nb_images_to_stack < LSF_LOW_STACK_WARN) {
		gchar *msg = g_strdup_printf(_("Only %d frames are selected, fewer than %d. The scatter of each pixel estimated "
					"from so few frames is unreliable, and any rejection decided on it stays in the stack."),
				stack_gui_args()->nb_images_to_stack, LSF_LOW_STACK_WARN);
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
	/* the floor and range frozen with the first batch (Section 10.3) */
	if (weighting == WFWHM_WEIGHT || weighting == NBSTARS_WEIGHT) {
		s->weight_floor = gtk_spin_button_get_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder,
						"stack_lsf_weight_floor")));
		s->weight_range = gtk_spin_button_get_value(GTK_SPIN_BUTTON(gtk_builder_get_object(gui.builder,
						"stack_lsf_weight_range")));
	}
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
	/* feathering, and the levels of the frames matched on their overlaps */
	a->feather = blend_feather();
	a->overlap_norm = blend_overlap_norm();
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


/* the Stacking button in a Living Stack mode */
void lsf_stacking_start(void) {
	if (stack_lsf_append_mode())
		start_lsf_append();
	else start_lsf_stacking();
}

/* The widgets of the tab for the result mode, from update_stack_interface()
 * once it set them for the stacking method. A Living Stack is a mean with
 * rejection made by lsf_init, which applies the registration itself with
 * its own framing, always writes float32, and has its own rejection maps
 * option (stack_lsf_rejmaps), as its maps are counts stored in the file,
 * not images. Returns TRUE when it changed the method, whose handler comes
 * back to update_stack_interface(), which then returns. */
gboolean lsf_stacking_update_widgets(int stack_method, gboolean can_reframe) {
	static GtkWidget *go_stack = NULL, *widgetnormalize = NULL, *force_norm = NULL, *output_norm = NULL,
			*RGB_equal = NULL, *fast_norm = NULL, *max_framing = NULL, *upscale_at_stacking = NULL,
			*blend_frame = NULL, *lsf_grid = NULL, *lsf_framing = NULL, *force32b = NULL, *rejmaps = NULL,
			*merge_rejmaps = NULL, *append_box = NULL, *output_box = NULL, *lsf_framing_label = NULL,
			*lsf_comp_label = NULL, *lsf_comp = NULL, *lsf_hist_same = NULL, *lsf_rejmaps = NULL,
			*weighing_combo = NULL, *siglow = NULL, *sighigh = NULL, *rejec_combo = NULL, *overlap_norm = NULL,
			*lsf_overlap_norm = NULL;
	static GtkDropDown *method_combo = NULL;
	if (!go_stack) {
		overlap_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_norm_overlap"));
		lsf_overlap_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_overlap_norm"));
		GtkWidget *feather = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_feather"));
		GtkWidget *feather_px = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_feather_px"));
		g_signal_connect(feather, "value-changed", G_CALLBACK(on_lsf_feather_changed), feather_px);
		on_lsf_feather_changed(GTK_SPIN_BUTTON(feather), feather_px);
		go_stack = GTK_WIDGET(gtk_builder_get_object(gui.builder, "gostack_button"));
		method_combo = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "comboboxstack_methods"));
		widgetnormalize = GTK_WIDGET(gtk_builder_get_object(gui.builder, "combonormalize"));
		force_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkforcenorm"));
		fast_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "checkfastnorm"));
		output_norm = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_normalise_to_max"));
		RGB_equal = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_RGBequal"));
		max_framing = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_maximize_framing"));
		upscale_at_stacking = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_upscale_at_stacking"));
		blend_frame = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_blend_frame"));
		lsf_grid = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_grid"));
		lsf_framing = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_framing"));
		force32b = GTK_WIDGET(gtk_builder_get_object(gui.builder, "check_force32b"));
		rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "rejmaps_checkbutton"));
		merge_rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "merge_rejmaps_checkbutton"));
		rejec_combo = GTK_WIDGET(gtk_builder_get_object(gui.builder, "comborejection"));
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
	}
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
	/* facts of the file are locked; the rejection may be overridden */
	GtkWidget *locked[] = { RGB_equal, output_norm, weighing_combo };
	for (int i = 0; i < (int) G_N_ELEMENTS(locked); i++)
		gtk_widget_set_sensitive(locked[i], !append);
	gtk_widget_set_sensitive(rejec_combo, TRUE);
	gtk_widget_set_sensitive(siglow, TRUE);
	gtk_widget_set_sensitive(sighigh, TRUE);
	if (append) {
		gtk_widget_set_sensitive(widgetnormalize, FALSE);
		gtk_widget_set_sensitive(fast_norm, FALSE);
	}
	if (!lsf) {
		/* back from a Living Stack mode: as on_comborejection_changed shows them */
		gboolean has_rej = gtk_drop_down_get_selected(GTK_DROP_DOWN(rejec_combo)) != NO_REJEC;
		gtk_widget_set_visible(rejmaps, has_rej);
		gtk_widget_set_visible(merge_rejmaps, has_rej);
		gtk_widget_set_visible(overlap_norm, TRUE);
		return FALSE;
	}
	if (stack_method != STACK_MEAN) {
		gtk_drop_down_set_selected(method_combo, STACK_MEAN);
		return TRUE;
	}
	gtk_widget_set_visible(max_framing, FALSE);
	gtk_widget_set_visible(upscale_at_stacking, FALSE);
	/* blending is in the Living Stack options, whatever the framing:
	 * feathering, and the normalization on the overlaps, with the Stack for
	 * an append, between the frames for a first batch of several panels */
	gtk_widget_set_visible(blend_frame, FALSE);
	gtk_widget_set_tooltip_text(lsf_overlap_norm, append ?
			_("For a mosaic panel: the frame of the batch that overlaps the Living Stack most is matched to it where "
				"they overlap, and the other frames to that one, instead of the normalization frozen by the first "
				"batch. The batch should show one panel. Only the level is matched, the scale stays the frozen one.") :
			_("For a first batch of several mosaic panels, plate-solved together: after their normalization by "
				"their global statistics, the levels of the frames are matched where they overlap, one frame after "
				"the other from the reference frame. Only the level is matched, the scale stays the frame's own."));
	gtk_widget_set_visible(rejmaps, FALSE);
	gtk_widget_set_visible(merge_rejmaps, FALSE);
	/* with the registration already applied, the frames define the grid */
	gtk_widget_set_sensitive(lsf_framing, can_reframe);
	lsf_stacking_update_hist_widgets();
	lsf_stacking_update_weight_widgets();
	/* an append can only add rejection maps to a file that has none */
	gtk_widget_set_visible(lsf_rejmaps, !append || (append_info && !append_info->file.settings.rejmaps));
	gtk_check_button_set_label(GTK_CHECK_BUTTON(lsf_rejmaps), append ? _("Add rejection maps to the file, counted from this append") :
			_("Keep rejection maps in the file"));
	return FALSE;
}

/* append_problems() and append_registration_note() read files: the frame
 * metadata, the file systems and the registration binding. They are kept
 * here and computed again only when what they depend on changes, as
 * update_stack_interface() runs on most widget changes. start_lsf_append()
 * reads the file again, so the checks it makes are always fresh. */
static struct {
	gchar *key;
	gchar *problems, *note;
	gboolean fatal, reusable;
} append_checks = { 0 };

void lsf_stacking_invalidate_checks(void) {
	g_clear_pointer(&append_checks.key, g_free);
}

static void update_append_checks(guint method) {
	GtkWidget *rejmaps = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_rejmaps"));
	gchar *out = append_out_path();
	GtkDropDown *disto = GTK_DROP_DOWN(gtk_builder_get_object(gui.builder, "stack_lsf_append_disto"));
	gchar *key = g_strdup_printf("%u|%s|%s|%d|%d|%d|%d|%u|%u|%d|%s", append_serial, com.wd ? com.wd : "",
			com.seq.seqname ? com.seq.seqname : "", com.seq.number, com.seq.selnum, com.seq.reference_image,
			stack_gui_args()->nb_images_to_stack, method, gtk_drop_down_get_selected(disto),
			gtk_widget_get_visible(rejmaps) && siril_toggle_get_active(rejmaps), out ? out : "");
	g_free(out);
	if (!g_strcmp0(key, append_checks.key)) {
		g_free(key);
		return;
	}
	g_free(append_checks.key);
	g_free(append_checks.problems);
	g_free(append_checks.note);
	append_checks.key = key;
	append_checks.problems = append_problems(&append_checks.fatal);
	append_checks.reusable = FALSE;
	/* a registration made beforehand is used by the star method only */
	append_checks.note = method == APPEND_REG_STARS ? append_registration_note(&append_checks.reusable) : NULL;
}

/* The append box, the frame count and the Stacking button in append mode,
 * from update_stack_interface() once it counted the frames; listed when the
 * stacking filters leave frames to append */
void lsf_stacking_update_append(gboolean listed) {
	static GtkWidget *go_stack = NULL, *append_framing = NULL, *expand_step = NULL, *expand_force = NULL,
			*info_label = NULL, *force_reg = NULL, *disto_label = NULL, *disto = NULL, *weight_inverse = NULL;
	static GtkLabel *result_label = NULL;
	if (!go_stack) {
		go_stack = GTK_WIDGET(gtk_builder_get_object(gui.builder, "gostack_button"));
		result_label = GTK_LABEL(gtk_builder_get_object(gui.builder, "stackfilter_label"));
		append_framing = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_framing"));
		expand_step = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_expand_step"));
		expand_force = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_expand_force"));
		info_label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_info_label"));
		force_reg = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_force_reg"));
		disto_label = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_disto_label"));
		disto = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_append_disto"));
		weight_inverse = GTK_WIDGET(gtk_builder_get_object(gui.builder, "stack_lsf_weight_inverse"));
	}
	gboolean grow = gtk_drop_down_get_selected(GTK_DROP_DOWN(append_framing)) == 0;
	guint method = append_regmethod();
	gboolean noreg = method == APPEND_REG_NONE;
	gtk_widget_set_sensitive(append_framing, !noreg);
	gtk_widget_set_sensitive(expand_step, grow && !noreg);
	gtk_widget_set_sensitive(expand_force, grow && !noreg);
	/* the distortion model, for an undistorted file only (Section 10.4) */
	gboolean undistorted = append_info && append_info->file.settings.undistorted;
	gtk_widget_set_visible(disto_label, undistorted);
	gtk_widget_set_visible(disto, undistorted);
	gtk_widget_set_sensitive(disto, !noreg);
	update_append_checks(method);
	gtk_widget_set_visible(force_reg, append_checks.reusable);
	/* the switch to inverse weighting, for a file whose weighting still has
	 * the range of its first batch (Section 10.3) */
	const lsf_settings *ws = append_info ? &append_info->file.settings : NULL;
	gboolean can_switch = ws && (ws->weighting == WFWHM_WEIGHT || ws->weighting == NBSTARS_WEIGHT)
		&& ws->weight_model == LSF_WMODEL_RANGE;
	gtk_widget_set_visible(weight_inverse, can_switch);
	GString *text = g_string_new(append_summary ? append_summary : _("Choose a Living Stack file."));
	if (append_checks.problems)
		g_string_append_printf(text, "\n%s", append_checks.problems);
	if (append_checks.note)
		g_string_append_printf(text, "\n%s", append_checks.note);
	if (can_switch && siril_toggle_get_active(weight_inverse))
		g_string_append(text, _("\nThis append switches the file to inverse weighting, for good: the frames already "
					"merged keep their weights."));
	gtk_label_set_text(GTK_LABEL(info_label), text->str);
	g_string_free(text, TRUE);
	gchar *labelbuffer = g_strdup_printf(_("Appending %d images of the %d of the sequence"),
			stack_gui_args()->nb_images_to_stack, com.seq.number);
	gtk_label_set_text(result_label, labelbuffer);
	g_free(labelbuffer);
	gtk_widget_set_sensitive(go_stack, append_info && !append_checks.fatal && listed);
}
