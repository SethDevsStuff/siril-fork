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

#include <gtk/gtk.h>
#include <string.h>
#include "core/siril.h"
#include "core/siril_log.h"
#include "core/processing.h"
#include "algos/siril_wcs.h"
#include "algos/photometry.h"
#include "algos/occult_time.h"
#include "algos/comparison_stars.h"
#include "gui-gtk4/occultation.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/utils.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/plot.h"
#include "gui-gtk4/image_display.h"
#include "io/sequence.h"

static GtkWidget *dialog = NULL;	// the window, a GtkWindow
static GtkWidget *delay_cam = NULL;
static GtkWidget *apply_offset = NULL;
static GtkWidget *res_data_grp, *pulse_lenght_grp, *data_grp;
static GtkWidget *ms_unit, *min_unit;
static GtkWidget *pls_len;
static double delay = 0.0;

static GtkWidget *occult_ok_button = NULL;

static void on_occult_response(GtkWindow *self, gint response_id, gpointer user_data);
static void on_find_clicked(GtkButton *button, gpointer user_data);

/* GtkDialog → GtkWindow: bridge per-button "clicked" signals back to the
 * legacy on_occult_response shape (cf. compstars.c). */
static void occult_btn_ok_clicked(GtkButton *btn, gpointer ud)    { (void)btn; on_occult_response(GTK_WINDOW(dialog), GTK_RESPONSE_ACCEPT, ud); }
static void occult_btn_close_clicked(GtkButton *btn, gpointer ud) { (void)btn; on_occult_response(GTK_WINDOW(dialog), GTK_RESPONSE_REJECT, ud); }

static void build_the_dialog() {
	/* GtkDialog → GtkWindow (GTK4). */
	dialog = gtk_window_new();
	gtk_window_set_title(GTK_WINDOW(dialog), _("Calibrate Timestamps over PPS"));
	gtk_window_set_default_size(GTK_WINDOW(dialog), 400, 400);
	gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
	gtk_window_set_modal(GTK_WINDOW(dialog), TRUE);
	/* A modal GtkWindow with no transient parent grabs input with no anchor
	 * and can wedge the whole session on Wayland/Xorg. */
	gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(lookup_widget("control_window")));
	gtk_window_set_hide_on_close(GTK_WINDOW(dialog), TRUE);
	g_signal_connect(G_OBJECT(dialog), "close-request", G_CALLBACK(siril_widget_hide_on_delete), NULL);

	// Definition of all the graphical items
	// Top label
	GtkWidget *label = gtk_label_new(_("Time re-Calibration of dataset a 1PPS signal trigger"));
	gtk_label_set_wrap(GTK_LABEL(label), TRUE);
	gtk_widget_set_margin_start(label, 15);
	gtk_widget_set_margin_end(label, 15);
	gtk_widget_set_margin_top(label, 15);
	gtk_widget_set_margin_bottom(label, 15);

	// Explanation text
	GtkWidget *label2 = gtk_label_new(_("The first step aims at finding the necessary delay to re-time\nand synchronize the timestamps to a 1PPS signal from a GPS module.\n\nTo get this parameter, you need:\n -a loaded sequence (from a .ser file)\n -this sequence must have a pseudo-star blinking at a beat of 1 pulse per second\n -make a selection around this pseudo-star, including enough background\n\nThen you can either click on Find delay or write your own delay if you know it\n"));
	gtk_label_set_wrap(GTK_LABEL(label2), TRUE);
	gtk_widget_set_margin_start(label2, 15);
	gtk_widget_set_margin_top(label2, 15);
	gtk_widget_set_margin_end(label2, 15);

	// The relevant data
	// Defines the results box, 3 lines
	res_data_grp = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
	gtk_box_set_homogeneous(GTK_BOX(res_data_grp), TRUE);
	gtk_widget_set_tooltip_text(res_data_grp, _("Find the offset delay in the loaded sequence if it embedds a 1-PPS blinking star"));

// Line#1: Lenght of the optical pulse
	// Definition of the 3 elements horizontal sub-box
	pulse_lenght_grp = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 3);
	gtk_box_set_homogeneous(GTK_BOX(pulse_lenght_grp), TRUE);
	gtk_widget_set_tooltip_text(pulse_lenght_grp, _("Used pulse width of the PPS"));
	// 1st element, a label
	GtkWidget *label1_pls_len = gtk_label_new(_("Pulse lenght: "));
	gtk_widget_set_halign(label1_pls_len, GTK_ALIGN_END);
	gtk_box_append(GTK_BOX(pulse_lenght_grp), label1_pls_len);

	// 2nd element, the pulse width entry
	pls_len = gtk_entry_new();
	gchar *pls_str = g_strdup_printf("%1.2lf", com.pref.phot_set.pulse_width);
	gtk_editable_set_text(GTK_EDITABLE(pls_len), pls_str);
	g_free(pls_str);
	gtk_entry_set_alignment(GTK_ENTRY(pls_len), 0.5);
	gtk_widget_set_tooltip_text(pls_len, _("Used optical pulse width. The value must match the one used during the acquisition."));
	gtk_widget_set_halign(pls_len, GTK_ALIGN_CENTER);
	gtk_widget_set_margin_top(pls_len, 0);
	gtk_widget_set_margin_bottom(pls_len, 0);
	gtk_box_append(GTK_BOX(pulse_lenght_grp), pls_len);
	// 3rd element, another label
	GtkWidget *label2_pls_unit = gtk_label_new(_("(ms)"));
	gtk_widget_set_halign(label2_pls_unit, GTK_ALIGN_START);
	gtk_box_append(GTK_BOX(pulse_lenght_grp), label2_pls_unit);
	// and finally include that box to the upper level box
	gtk_box_append(GTK_BOX(res_data_grp), pulse_lenght_grp);

// Line#2: The result fields
	//// A line with 3 elements
	data_grp = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 3);
	gtk_box_set_homogeneous(GTK_BOX(data_grp), TRUE);
	gtk_widget_set_tooltip_text(data_grp, _("Output data of the time synchronization process or user setting for the next step"));
	//// Item #1, the label
	GtkWidget *label_delay_cam = gtk_label_new(_("Computed Delay:"));
	gtk_widget_set_halign(label_delay_cam, GTK_ALIGN_END);
	gtk_widget_set_margin_start(label_delay_cam, 15);
	gtk_box_append(GTK_BOX(data_grp), label_delay_cam);
	//// Item#2: the Delay field
	delay_cam = gtk_entry_new();
	gtk_editable_set_text(GTK_EDITABLE(delay_cam), "0.0");
	gtk_entry_set_alignment(GTK_ENTRY(delay_cam), 0.5);
	gtk_widget_set_tooltip_text(delay_cam, _("Delay"));
	gtk_widget_set_halign(delay_cam, GTK_ALIGN_CENTER);
	gtk_box_append(GTK_BOX(data_grp), delay_cam);
	//// Item#3: the unit button, ms or min
	GtkWidget *unit_box;
	unit_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
	gtk_box_set_homogeneous(GTK_BOX(unit_box), TRUE);
	gtk_widget_set_tooltip_text(unit_box, _("(ms) is usefull for occultation resync whereas (min) is usefull for transit resync"));
	/* GtkRadioButton removed in GTK4; use GtkCheckButton + group. */
	ms_unit = gtk_check_button_new_with_label(_("ms"));
	min_unit = gtk_check_button_new_with_label(_("min"));
	gtk_check_button_set_group(GTK_CHECK_BUTTON(min_unit), GTK_CHECK_BUTTON(ms_unit));
	gtk_check_button_set_active(GTK_CHECK_BUTTON(ms_unit), TRUE);
	gtk_box_append(GTK_BOX(unit_box), ms_unit);
	gtk_box_append(GTK_BOX(unit_box), min_unit);
	gtk_box_append(GTK_BOX(data_grp), unit_box);
	gtk_box_append(GTK_BOX(res_data_grp), data_grp);	// Final Add to the upper container

// Line#3: The button
	GtkWidget *find_delay = gtk_button_new_with_mnemonic(_("Find Delay"));
	g_signal_connect(find_delay, "clicked", G_CALLBACK(on_find_clicked), NULL);
	gtk_box_append(GTK_BOX(res_data_grp), find_delay);	// Final Add to the upper container

	GtkWidget *sep = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
	gtk_widget_set_size_request(sep, 1, 5);
	gtk_widget_set_margin_top(sep, 15);

	GtkWidget *label3 = gtk_label_new(_("Once the delay is set, you can choose to apply it or not \n"));
	gtk_widget_set_halign(label3, GTK_ALIGN_START);
	gtk_widget_set_margin_start(label3, 15);
	gtk_widget_set_margin_end(label3, 15);
	gtk_widget_set_margin_top(label3, 15);

	apply_offset = gtk_check_button_new_with_label(_("Apply Offset"));
	gtk_widget_set_halign(apply_offset, GTK_ALIGN_CENTER);
	gtk_widget_set_margin_bottom(apply_offset, 35);
	gtk_widget_set_margin_top(apply_offset, 5);
	gtk_widget_set_tooltip_text(apply_offset, _("The previously processed or set value can be applied the the loaded sequence.\nHence, the new timestamps shall be in the form:\n	new_timestamp = old_timestamp + offset delay"));

// Gather the graphic items
	GtkWidget *content_area = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
	gtk_box_append(GTK_BOX(content_area), label);
	gtk_box_append(GTK_BOX(content_area), label2);
	gtk_box_append(GTK_BOX(content_area), res_data_grp);
	gtk_box_append(GTK_BOX(content_area), sep);
	gtk_box_append(GTK_BOX(content_area), label3);
	gtk_box_append(GTK_BOX(content_area), apply_offset);

	/* Action area: OK / Close buttons (GtkWindow has no built-in action area). */
	GtkWidget *bbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
	gtk_widget_set_halign(bbox, GTK_ALIGN_END);
	gtk_widget_set_margin_end(bbox, 12);
	gtk_widget_set_margin_bottom(bbox, 12);
	GtkWidget *btn_close = gtk_button_new_with_mnemonic(_("_Close"));
	g_signal_connect(btn_close, "clicked", G_CALLBACK(occult_btn_close_clicked), NULL);
	gtk_box_append(GTK_BOX(bbox), btn_close);
	occult_ok_button = gtk_button_new_with_mnemonic(_("_OK"));
	gtk_widget_add_css_class(occult_ok_button, "suggested-action");
	g_signal_connect(occult_ok_button, "clicked", G_CALLBACK(occult_btn_ok_clicked), NULL);
	gtk_box_append(GTK_BOX(bbox), occult_ok_button);
	gtk_box_append(GTK_BOX(content_area), bbox);
	gtk_window_set_default_widget(GTK_WINDOW(dialog), occult_ok_button);

	gtk_window_set_child(GTK_WINDOW(dialog), content_area);
}

// the public getter
GtkWidget *get_occult_dialog() {
	if (!dialog)
		build_the_dialog();

	return dialog;
}

/* The processing thread.  Lives in the GUI module (not in libsiril/occult_time.c)
 * because it drives GUI refreshes (queue_redraw) and schedules the GUI idle
 * callback; the pure timing math stays in algos/occult_time.c. */
gpointer occultation_worker(gpointer arg) {
	int retval = 0;
	gboolean tmp_t_delayed = com.pref.phot_set.t_delayed;
	struct light_curve_args *args = (struct light_curve_args *)arg;
	framing_mode framing = REGISTERED_FRAME;

	siril_log_info(_("Entering the Occulation Time precedure.\n"));
	siril_log_info(_("PSF of blinking-star...\n"));

	com.pref.phot_set.t_delayed = TRUE;	// workaround to avoid all the warning/error messages in seqpsf...
	if (seqpsf(args->seq, args->layer, FALSE, TRUE, FALSE, framing, FALSE, TRUE)) {
		siril_log_error(_("Something went wrong with the PSF analisys. Try to enlarge the selection or check the blinking star is not too faint\n"));
		return GINT_TO_POINTER(retval);
	}
	com.pref.phot_set.t_delayed = tmp_t_delayed;	// ...end of the workaround

	if (args->seq == &com.seq)
		queue_redraw(REDRAW_OVERLAY);

	memset(&com.selection, 0, sizeof(rectangle));

	/* analyse data and compute the offset */
	if (!retval)
		retval = time_resync(args);

	siril_add_idle(end_occultation_worker, args);
	return GINT_TO_POINTER(retval);
}

static void on_find_clicked(GtkButton *button, gpointer user_data) {
	(void)button; (void)user_data;
	if (com.seq.photometry[0] != NULL) free_photometry_set(&com.seq, 0);
	control_window_switch_to_tab(OUTPUT_LOGS);

	int layer = -1;
	if (com.seq.regparam) {
		for (int i = 0; i < com.seq.nb_layers; i++) {
			if (com.seq.regparam[i]) {
				layer = i;
				break;
			}
		}
	}
	if (layer == -1) {
		siril_log_debug("unregistered sequence\n");
		//siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("The sequence must be registered"));
		if (com.seq.nb_layers == 3)
			layer = 1;
		else layer = 0;
	}

	// An area must be drawn around the blinking star
	if (com.selection.w == 0 || com.selection.h  == 0) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("Error"), _("Make a selection aroud the blinking star first"));
		return;
	}
	// This area must be large enought wrt the current photometric aperture (outer ring)
	if ((double)com.selection.w <= com.pref.phot_set.outer || (double)com.selection.h <= com.pref.phot_set.outer) {
		siril_log_warning(_("The selection has benen resized \n"));
		com.selection.w = (int)(1.1 * com.pref.phot_set.outer);
		com.selection.h = (int)(1.1 * com.pref.phot_set.outer);
	}

	com.pref.phot_set.pulse_width = g_strtod(gtk_editable_get_text(GTK_EDITABLE(pls_len)), NULL);	// Use the entered pulse width by the user if it's been modified

	struct light_curve_args *args = calloc(1, sizeof(struct light_curve_args));
	args->seq = &com.seq;
	args->layer = layer;

	start_in_new_thread(occultation_worker, args);

	delay = 0.0;
	gchar *delay_str = g_strdup_printf("%0.3lf", delay);
	gtk_editable_set_text(GTK_EDITABLE(delay_cam), delay_str);
	g_free(delay_str);
	gtk_check_button_set_active(GTK_CHECK_BUTTON(min_unit), FALSE);
	gtk_check_button_set_active(GTK_CHECK_BUTTON(ms_unit), TRUE);
	control_window_switch_to_tab(OUTPUT_LOGS);
}

gboolean end_occultation_worker(gpointer p) {
	if (!com.script) {
		struct light_curve_args *args = (struct light_curve_args *)p;
		args->seq = &com.seq;
		delay = com.pref.phot_set.time_offset;
		gchar *delay_str = g_strdup_printf("%0.3lf", delay);
		gtk_editable_set_text(GTK_EDITABLE(delay_cam), delay_str);
		g_free(delay_str);
		control_window_switch_to_tab(OUTPUT_LOGS);
		free_light_curve_args(args);
	}
	if (sequence_is_loaded()) {
		drawPlot();
		notify_new_photometry();	// switch to and update plot tab
		redraw(REDRAW_OVERLAY);
	}
	return end_generic(NULL);
}

static void on_occult_response(GtkWindow *self, gint response_id, gpointer user_data) {
	(void)self; (void)user_data;
	siril_log_debug("got response event\n");
	if (response_id != GTK_RESPONSE_ACCEPT) {
		gtk_widget_set_visible(dialog, FALSE);
		reactivate_parent(dialog);
		return;
	}

	control_window_switch_to_tab(OUTPUT_LOGS);
	gboolean use_offset = siril_toggle_get_active(GTK_WIDGET(apply_offset));
	gboolean use_ms = siril_toggle_get_active(GTK_WIDGET(ms_unit));
	if (use_offset) {
		com.pref.phot_set.t_delayed = TRUE;
		double dbl_offset = g_strtod(gtk_editable_get_text(GTK_EDITABLE(delay_cam)), NULL);
		com.pref.phot_set.time_offset = use_ms ? dbl_offset : 60000.0 * dbl_offset;	// Converted in (ms)
		siril_log_message(_("Applied offset: %0.3lf (%s) \n"), dbl_offset, use_ms ? "ms" : "min");
	}
	else {
		siril_log_message(_("No offset applied \n"));
		com.pref.phot_set.t_delayed = FALSE;
		com.pref.phot_set.time_offset = 0.0;
	}
	gtk_widget_set_visible(dialog, FALSE);
	reactivate_parent(dialog);
}
