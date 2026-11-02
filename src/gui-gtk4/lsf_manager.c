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

/* The Living Stack dialog: one LSF file's summary, ledger and
 * history, and the lsf_util, lsf_export and lsf_maps jobs on it. Everything
 * that changes the file goes through the same jobs as the commands. */

#include <math.h>
#include <string.h>
#include <glib/gstdio.h>

#include "algos/siril_wcs.h"
#include "core/processing.h"
#include "core/proto.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "gui-gtk4/dialogs.h"
#include "gui-gtk4/lsf_manager.h"
#include "gui-gtk4/message_dialog.h"
#include "gui-gtk4/progress_and_log.h"
#include "gui-gtk4/utils.h"
#include "io/image_format_fits.h"
#include "io/single_image.h"
#include "lsf/lsf.h"
#include "lsf/lsf_commands.h"
#include "lsf/lsf_fits.h"
#include "lsf/lsf_settings.h"
#include "lsf/lsf_wcs.h"

static lsf_info *info = NULL;		// the file shown, or NULL
static gboolean stack_edited = FALSE;	// the last check of this file found check 8 failed
static gchar *pending_out = NULL;	// export or maps output of the running job, for its result

/* ── ledger rows ───────────────────────────────────────────────────────── */

#define SIRIL_TYPE_LSF_ROW (siril_lsf_row_get_type())
G_DECLARE_FINAL_TYPE(SirilLsfRow, siril_lsf_row, SIRIL, LSF_ROW, GObject)
struct _SirilLsfRow {
	GObject parent_instance;
	gint index, batch, nstars;
	gchar *date, *file, *filter;
	gdouble exptime, fwhm, wfwhm, weight, noise;
	gboolean at_floor;	// weight at LSF_WEIGHT_FLOOR (Section 10.3)
};
G_DEFINE_TYPE(SirilLsfRow, siril_lsf_row, G_TYPE_OBJECT)
static void siril_lsf_row_finalize(GObject *obj) {
	SirilLsfRow *r = SIRIL_LSF_ROW(obj);
	g_free(r->date);
	g_free(r->file);
	g_free(r->filter);
	G_OBJECT_CLASS(siril_lsf_row_parent_class)->finalize(obj);
}
static void siril_lsf_row_class_init(SirilLsfRowClass *klass) {
	G_OBJECT_CLASS(klass)->finalize = siril_lsf_row_finalize;
}
static void siril_lsf_row_init(SirilLsfRow *self) { (void)self; }

typedef enum { COL_INDEX, COL_BATCH, COL_DATE, COL_FILE, COL_FILTER, COL_EXPTIME, COL_FWHM, COL_WFWHM,
	COL_NSTARS, COL_WEIGHT, COL_NOISE } LsfColKind;

static GListStore *ledger_store = NULL;
static GtkColumnView *ledger_view = NULL;

static void setup_label_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	GtkWidget *lbl = gtk_label_new(NULL);
	gtk_label_set_xalign(GTK_LABEL(lbl), GPOINTER_TO_INT(u) ? 1.0 : 0.0);
	gtk_list_item_set_child(li, lbl);
}

static void bind_cb(GtkSignalListItemFactory *f, GtkListItem *li, gpointer u) {
	LsfColKind kind = (LsfColKind) GPOINTER_TO_INT(u);
	GtkLabel *lbl = GTK_LABEL(gtk_list_item_get_child(li));
	SirilLsfRow *r = SIRIL_LSF_ROW(gtk_list_item_get_item(li));
	if (!r)
		return;
	gchar *txt = NULL;
	switch (kind) {
	case COL_INDEX: txt = g_strdup_printf("%d", r->index); break;
	case COL_BATCH: txt = g_strdup_printf("%d", r->batch); break;
	case COL_DATE: txt = g_strdup(r->date ? r->date : ""); break;
	case COL_FILE: txt = g_strdup(r->file ? r->file : ""); break;
	case COL_FILTER: txt = g_strdup(r->filter ? r->filter : ""); break;
	case COL_EXPTIME: txt = g_strdup_printf("%g", r->exptime); break;
	case COL_FWHM: txt = r->fwhm > 0.0 ? g_strdup_printf("%.2f", r->fwhm) : g_strdup(""); break;
	case COL_WFWHM: txt = r->wfwhm > 0.0 ? g_strdup_printf("%.2f", r->wfwhm) : g_strdup(""); break;
	case COL_NSTARS: txt = g_strdup_printf("%d", r->nstars); break;
	case COL_WEIGHT: txt = r->at_floor ? g_strdup_printf(_("%.4f (floor)"), r->weight) : g_strdup_printf("%.4f", r->weight); break;
	case COL_NOISE: txt = g_strdup_printf("%.3e", r->noise); break;
	}
	gtk_label_set_text(lbl, txt);
	g_free(txt);
}

#define DEFINE_CMP(name, expr) \
	static gint cmp_##name(gconstpointer a, gconstpointer b, gpointer u) { \
		const SirilLsfRow *x = a, *y = b; (void)u; return expr; }
#define NUM_CMP(f) (x->f < y->f) ? -1 : (x->f > y->f)
DEFINE_CMP(index, NUM_CMP(index))
DEFINE_CMP(batch, NUM_CMP(batch))
DEFINE_CMP(date, g_strcmp0(x->date, y->date))
DEFINE_CMP(file, g_strcmp0(x->file, y->file))
DEFINE_CMP(filter, g_strcmp0(x->filter, y->filter))
DEFINE_CMP(exptime, NUM_CMP(exptime))
DEFINE_CMP(fwhm, NUM_CMP(fwhm))
DEFINE_CMP(wfwhm, NUM_CMP(wfwhm))
DEFINE_CMP(nstars, NUM_CMP(nstars))
DEFINE_CMP(weight, NUM_CMP(weight))
DEFINE_CMP(noise, NUM_CMP(noise))
#undef NUM_CMP
#undef DEFINE_CMP

static void add_column(const char *title, LsfColKind kind, gboolean numeric, GCompareDataFunc cmp) {
	GtkSignalListItemFactory *f = GTK_SIGNAL_LIST_ITEM_FACTORY(gtk_signal_list_item_factory_new());
	g_signal_connect(f, "setup", G_CALLBACK(setup_label_cb), GINT_TO_POINTER(numeric));
	g_signal_connect(f, "bind", G_CALLBACK(bind_cb), GINT_TO_POINTER(kind));
	GtkColumnViewColumn *c = gtk_column_view_column_new(title, GTK_LIST_ITEM_FACTORY(f));
	gtk_column_view_column_set_resizable(c, TRUE);
	gtk_column_view_column_set_sorter(c, GTK_SORTER(gtk_custom_sorter_new(cmp, NULL, NULL)));
	gtk_column_view_append_column(ledger_view, c);
	g_object_unref(c);
}

static void build_ledger_view() {
	if (ledger_view)
		return;
	ledger_store = g_list_store_new(SIRIL_TYPE_LSF_ROW);
	ledger_view = GTK_COLUMN_VIEW(gtk_column_view_new(NULL));
	gtk_widget_add_css_class(GTK_WIDGET(ledger_view), "siril-dense-rows");
	add_column("#", COL_INDEX, TRUE, cmp_index);
	add_column(_("Batch"), COL_BATCH, TRUE, cmp_batch);
	add_column(_("Date"), COL_DATE, FALSE, cmp_date);
	add_column(_("File"), COL_FILE, FALSE, cmp_file);
	add_column(_("Filter"), COL_FILTER, FALSE, cmp_filter);
	add_column(_("Exp. (s)"), COL_EXPTIME, TRUE, cmp_exptime);
	add_column("FWHM", COL_FWHM, TRUE, cmp_fwhm);
	add_column("wFWHM", COL_WFWHM, TRUE, cmp_wfwhm);
	add_column(_("Stars"), COL_NSTARS, TRUE, cmp_nstars);
	add_column(_("Weight"), COL_WEIGHT, TRUE, cmp_weight);
	add_column(_("Noise"), COL_NOISE, TRUE, cmp_noise);
	GtkSorter *sorter = gtk_column_view_get_sorter(ledger_view);
	GtkSortListModel *sorted = gtk_sort_list_model_new(G_LIST_MODEL(g_object_ref(ledger_store)),
			sorter ? g_object_ref(sorter) : NULL);
	gtk_column_view_set_model(ledger_view, GTK_SELECTION_MODEL(gtk_no_selection_new(G_LIST_MODEL(sorted))));
	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(gtk_builder_get_object(gui.builder, "lsfm_ledger_scrolled")),
			GTK_WIDGET(ledger_view));
}

/* ── summary ───────────────────────────────────────────────────────────── */

static gchar *short_date(const gchar *iso) {
	return iso ? g_strndup(iso, 10) : g_strdup("?");
}

gchar *lsf_gui_describe(const lsf_info *in) {
	const lsf_file *f = &in->file;
	const lsf_settings *s = &f->settings;
	gchar *first = short_date(in->first_date), *last = short_date(in->last_date),
	      *added = short_date(in->last_added), *size = g_format_size(in->file_bytes);
	gchar *grid = s->oversample > 1 ? g_strdup_printf(_("%dx%d grid, drizzled x%d"), f->rx, f->ry, s->oversample) :
		g_strdup_printf(_("%dx%d grid"), f->rx, f->ry);
	gchar *expmap = f->has_exp ? g_strdup(_("exposure map kept")) : isfinite(f->ledger_exptime) ?
		g_strdup_printf(_("every frame %g s"), f->ledger_exptime) : g_strdup(_("exposure per pixel unknown"));
	gchar *msg = g_strdup_printf(_("%ld frames in %d batches, %.1f h of exposure (%s), frames from %s to %s, last "
				"appended %s.\n%s, %s, filter %s. Compression %s, rejection maps %s, %s."),
			f->ledger_rows, f->nb_batches, in->exposure / 3600.0, expmap, first, last, added, grid,
			f->nchan == 3 ? _("colour") : _("mono"), s->filter && s->filter[0] ? s->filter : _("none"),
			lsf_compression_name(f->compression), s->rejmaps ? _("kept") : _("not kept"), size);
	g_free(expmap);
	g_free(first);
	g_free(last);
	g_free(added);
	g_free(size);
	g_free(grid);
	/* the WCS, which registration through the WCS needs, and the distortion
	 * correction of the grid (Section 10.4) */
	gchar *wcs = has_wcs((fits *) &f->primary) ? lsf_wcs_describe(f->primary.keywords.wcslib, f->rx, f->ry)
		: g_strdup(_("none"));
	gchar *more;
	if (s->undistorted) {
		gchar *model = lsf_disto_describe(&s->disto);
		more = g_strdup_printf(_("%s\nWCS: %s. Undistorted grid: every append undistorts its frames, with %s by "
					"default."), msg, wcs, model);
		g_free(model);
	} else more = g_strdup_printf(_("%s\nWCS: %s. Frames are not undistorted."), msg, wcs);
	g_free(wcs);
	g_free(msg);
	/* the weighting, and for wfwhm and nbstars how it scales the frames
	 * (Section 10.3) */
	gchar *weights;
	if (s->weighting == WFWHM_WEIGHT || s->weighting == NBSTARS_WEIGHT) {
		gchar *scale = s->weight_model == LSF_WMODEL_INVERSE ? g_strdup(_("inverse, no range"))
			: s->weighting == WFWHM_WEIGHT ? g_strdup_printf(_("range of wFWHM %.2f to %.2f"), s->fwhm_min, s->fwhm_max)
			: g_strdup_printf(_("range of %.0f to %.0f stars"), s->stars_min, s->stars_max);
		weights = g_strdup_printf(_("%s\nWeighting %s: %s, floor %g."), more, lsf_weight_name(s->weighting), scale,
				s->weight_floor);
		g_free(scale);
	} else weights = g_strdup_printf(_("%s\nWeighting %s."), more, lsf_weight_name(s->weighting));
	g_free(more);
	return weights;
}

static GtkWidget *w(const char *id) {
	return GTK_WIDGET(gtk_builder_get_object(gui.builder, id));
}

/* the weight a wfwhm or nbstars frame gets at the floor, after the frozen
 * normalization (Section 10.3), or -1 when the scheme has no floor */
static double floor_weight(const lsf_settings *s) {
	if (s->weighting == WFWHM_WEIGHT && s->fwhm_norm > 0.0)
		return s->weight_floor / s->fwhm_norm;
	if (s->weighting == NBSTARS_WEIGHT && s->stars_norm > 0.0)
		return s->weight_floor / s->stars_norm;
	return -1.0;
}

static void fill_ledger() {
	g_list_store_remove_all(ledger_store);
	GString *batches = g_string_new(NULL);
	GPtrArray *rows = info ? lsf_read_ledger_view(info) : NULL;
	if (rows) {
		const lsf_settings *s = &info->file.settings;
		double wfloor = floor_weight(s);
		int nb = info->file.nb_batches, nchan = info->file.nchan;
		for (int b = 1; b <= nb; b++) {
			int n = 0, at_floor = 0, nfwhm = 0;
			double exp = 0.0, fwhm = 0.0;
			for (guint i = 0; i < rows->len; i++) {
				lsf_ledger_view *v = g_ptr_array_index(rows, i);
				if (v->batch != b)
					continue;
				n++;
				exp += v->exptime;
				if (v->fwhm > 0.0f) {
					fwhm += v->fwhm;
					nfwhm++;
				}
				at_floor += wfloor > 0.0 && v->weight[0] <= wfloor * (1.0 + 1e-6);
			}
			g_string_append_printf(batches, _("Batch %d: %d frames, %.0f s"), b, n, exp);
			if (nfwhm)
				g_string_append_printf(batches, _(", mean FWHM %.2f"), fwhm / nfwhm);
			if (at_floor)
				g_string_append_printf(batches, _(", %d at the weight floor"), at_floor);
			g_string_append(batches, ".\n");
		}
		if (wfloor > 0.0)
			g_string_append(batches, _("A weight at the floor marks a frame as bad as, or worse than, the worst frame "
						"of the first batch."));
		for (guint i = 0; i < rows->len; i++) {
			lsf_ledger_view *v = g_ptr_array_index(rows, i);
			SirilLsfRow *r = g_object_new(SIRIL_TYPE_LSF_ROW, NULL);
			r->index = i + 1;
			r->batch = v->batch;
			r->date = v->date_obs ? g_strndup(v->date_obs, 19) : NULL;
			r->file = v->filename ? g_path_get_basename(v->filename) : NULL;
			r->filter = g_strdup(v->filter);
			r->exptime = v->exptime;
			r->fwhm = v->fwhm;
			r->wfwhm = v->wfwhm;
			r->nstars = v->nstars;
			/* one weight and noise for mono, the green ones for colour */
			r->weight = v->weight[nchan == 3 ? 1 : 0];
			r->noise = v->noise[nchan == 3 ? 1 : 0];
			r->at_floor = wfloor > 0.0 && v->weight[0] <= wfloor * (1.0 + 1e-6);
			g_list_store_append(ledger_store, r);
			g_object_unref(r);
		}
		g_ptr_array_free(rows, TRUE);
	}
	gtk_label_set_text(GTK_LABEL(w("lsfm_batches")), batches->str);
	g_string_free(batches, TRUE);
}

static void fill_history() {
	GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(w("lsfm_history")));
	GPtrArray *cards = info ? lsf_read_history(info->file.path) : NULL;
	GString *text = g_string_new(NULL);
	for (guint i = 0; cards && i < cards->len; i++) {
		g_string_append(text, g_ptr_array_index(cards, i));
		g_string_append_c(text, '\n');
	}
	if (cards)
		g_ptr_array_free(cards, TRUE);
	gtk_text_buffer_set_text(buf, text->len ? text->str : _("No history recorded."), -1);
	g_string_free(text, TRUE);
}

/* reads the file of the entry and shows it */
static void load_file() {
	if (info) {
		lsf_info_free(info);
		g_free(info);
		info = NULL;
	}
	gchar *path = g_strstrip(g_strdup(gtk_editable_get_text(GTK_EDITABLE(w("lsfm_entry")))));
	gchar *summary = NULL;
	if (!path[0])
		summary = g_strdup(_("Choose a Living Stack file."));
	else if (!g_file_test(path, G_FILE_TEST_IS_REGULAR))
		summary = g_strdup(_("This file does not exist."));
	else {
		lsf_info *in = g_new0(lsf_info, 1);
		if (lsf_read_info(path, in)) {
			summary = g_strdup_printf(_("This is not a Living Stack that can be used: %s."), in->problem);
			lsf_info_free(in);
			g_free(in);
		} else {
			info = in;
			summary = lsf_gui_describe(in);
			if (in->no_append) {
				gchar *more = g_strdup_printf(_("%s\nThis version cannot append to it: %s."), summary, in->no_append);
				g_free(summary);
				summary = more;
			}
		}
	}
	gtk_label_set_text(GTK_LABEL(w("lsfm_summary")), summary);
	g_free(summary);

	gboolean ok = info != NULL;
	const char *actions[] = { "lsfm_check", "lsfm_compress", "lsfm_compression", "lsfm_copy", "lsfm_export",
		"lsfm_import_wcs", "lsfm_maps_save" };
	for (int i = 0; i < (int) G_N_ELEMENTS(actions); i++)
		gtk_widget_set_sensitive(w(actions[i]), ok);
	gtk_widget_set_sensitive(w("lsfm_repair"), ok && stack_edited);
	gtk_widget_set_sensitive(w("lsfm_adopt"), ok && stack_edited);
	gboolean rej = ok && info->file.settings.rejmaps;
	gtk_widget_set_sensitive(w("lsfm_map_rejlo"), rej);
	gtk_widget_set_sensitive(w("lsfm_map_rejhi"), rej);
	/* exposure per pixel: the map, or COUNT times the one exposure (Section 10.8) */
	gtk_widget_set_sensitive(w("lsfm_map_exp"), ok && (info->file.has_exp || isfinite(info->file.ledger_exptime)));
	if (ok) {
		gtk_drop_down_set_selected(GTK_DROP_DOWN(w("lsfm_compression")), info->file.compression == LSF_COMP_GZIP2 ? 1 :
				info->file.compression == LSF_COMP_GZIP1 ? 2 : 0);
		GtkEditable *prefix = GTK_EDITABLE(w("lsfm_maps_prefix"));
		if (!gtk_editable_get_text(prefix)[0] && g_str_has_suffix(path, ".lsf.fits")) {
			gchar *p = g_strndup(path, strlen(path) - strlen(".lsf.fits"));
			gtk_editable_set_text(prefix, p);
			g_free(p);
		}
	}
	g_free(path);
	fill_ledger();
	fill_history();
}

void on_lsfm_entry_changed(GtkEditable *editable, gpointer user_data) {
	stack_edited = FALSE;	// a check is about one file
	gtk_label_set_text(GTK_LABEL(w("lsfm_result")), "");
	gtk_editable_set_text(GTK_EDITABLE(w("lsfm_maps_prefix")), "");
	load_file();
}

void on_lsfm_browse_clicked(GtkButton *button, gpointer user_data) {
	SirilFileChooser *fc = siril_fc_open(GTK_WINDOW(w("lsf_manager_dialog")), GTK_FILE_CHOOSER_ACTION_OPEN);
	siril_fc_set_current_folder_path(fc, com.wd);
	siril_fc_set_select_multiple(fc, FALSE);
	siril_fc_add_filter_pattern(fc, _("Living Stack files (*.lsf.fits)"), "*.lsf.fits;*.LSF.FITS", TRUE);
	if (siril_fc_run(fc) == GTK_RESPONSE_ACCEPT) {
		gchar *file = siril_fc_get_filename(fc);
		if (file) {
			gtk_editable_set_text(GTK_EDITABLE(w("lsfm_entry")), file);
			gtk_editable_set_position(GTK_EDITABLE(w("lsfm_entry")), -1);
			g_free(file);
		}
	}
	siril_fc_destroy(fc);
}

/* ── jobs ──────────────────────────────────────────────────────────────── */

/* a file left for the user to deal with when replacing it failed (Section
 * 10.10): the message of the job says which files and what to do */
void lsf_gui_report_write(const lsf_job_result *res) {
	if (res->write_status == LSF_WRITE_KEPT_TMP)
		siril_message_dialog(GTK_MESSAGE_WARNING, _("The Living Stack was not updated"),
				res->message ? res->message : _("The file could not be replaced, see the log."));
	else if (res->write_status == LSF_WRITE_LEFT_AS_BAK)
		siril_message_dialog(GTK_MESSAGE_ERROR, _("The Living Stack must be renamed back"),
				res->message ? res->message : _("The file could not be replaced, see the log."));
}

/* Repair, Adopt or Compression written to a new file: the dialog then shows
 * it. A copy leaves the dialog on its source. */
static gboolean show_new_file = FALSE;

static void end_util(const lsf_job_result *res) {
	set_cursor_waiting(FALSE);
	lsf_gui_report_write(res);
	if (res->mode == LSF_CMD_UTIL)
		stack_edited = res->stack_modified;
	gchar *text = res->message ? g_strdup(res->message) :
		g_strdup(res->retval ? _("The operation failed, see the log.") : _("Done."));
	if (!res->retval && res->out && show_new_file) {
		/* loaded by the changed handler, which clears the result */
		gchar *path = g_canonicalize_filename(res->out, NULL);
		gtk_editable_set_text(GTK_EDITABLE(w("lsfm_entry")), path);
		gtk_editable_set_position(GTK_EDITABLE(w("lsfm_entry")), -1);
		g_free(path);
	} else load_file();
	gtk_label_set_text(GTK_LABEL(w("lsfm_result")), text);
	g_free(text);
}

static void end_export(const lsf_job_result *res) {
	set_cursor_waiting(FALSE);
	gchar *text = res->retval ? g_strdup(_("The export failed, see the log.")) :
		g_strdup_printf(_("The stack was exported to %s."), pending_out);
	gtk_label_set_text(GTK_LABEL(w("lsfm_result")), text);
	g_free(text);
}

static const struct { const char *id, *suffix; guint flag; } map_kinds[] = {
	{ "lsfm_map_neff", "neff", LSF_MAP_NEFF }, { "lsfm_map_noise", "noise", LSF_MAP_NOISE },
	{ "lsfm_map_sigma", "sigma", LSF_MAP_SIGMA }, { "lsfm_map_snr", "snr", LSF_MAP_SNR },
	{ "lsfm_map_rejlo", "rejlo", LSF_MAP_REJLO }, { "lsfm_map_rejhi", "rejhi", LSF_MAP_REJHI },
	{ "lsfm_map_exp", "exp", LSF_MAP_EXP }
};

static void end_maps(const lsf_job_result *res) {
	set_cursor_waiting(FALSE);
	if (res->retval) {
		gtk_label_set_text(GTK_LABEL(w("lsfm_maps_result")), _("The maps were not saved, see the log."));
		return;
	}
	gchar *text = g_strdup_printf(_("Maps saved as %s_*."), pending_out);
	gtk_label_set_text(GTK_LABEL(w("lsfm_maps_result")), text);
	g_free(text);
	if (!siril_toggle_get_active(w("lsfm_maps_open")))
		return;
	for (int i = 0; i < (int) G_N_ELEMENTS(map_kinds); i++) {
		GtkWidget *check = w(map_kinds[i].id);
		if (gtk_widget_get_sensitive(check) && siril_toggle_get_active(check)) {
			gchar *base = g_strdup_printf("%s_%s", pending_out, map_kinds[i].suffix);
			gchar *name = set_right_extension(base);
			if (name)
				open_single_image(name);
			g_free(name);
			g_free(base);
			break;
		}
	}
}

static void start_job(struct lsf_cmd_args *a) {
	gchar *cmd = lsf_cmd_args_to_command(a);
	if (cmd) {
		siril_log_message(_("Living Stack, equivalent command: %s\n"), cmd);
		g_free(cmd);
	}
	set_cursor_waiting(TRUE);
	if (lsf_start(a)) {
		PRINT_ANOTHER_THREAD_RUNNING;
		set_cursor_waiting(FALSE);
	}
}

/* an lsf_util job on the file shown */
static struct lsf_cmd_args *new_util_args() {
	if (!info || processing_is_job_active()) {
		if (info)
			PRINT_ANOTHER_THREAD_RUNNING;
		return NULL;
	}
	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_UTIL);
	if (!a)
		return NULL;
	a->path = g_strdup(info->file.path);
	a->end_func = end_util;
	show_new_file = FALSE;
	gtk_label_set_text(GTK_LABEL(w("lsfm_result")), _("Working…"));
	return a;
}

static gchar *choose_save_name(const char *title_filter, const char *pattern, const char *suggested);

/* the file name without .lsf.fits */
static gchar *lsf_stem(const char *path) {
	gchar *base = g_path_get_basename(path);
	gchar *stem = g_str_has_suffix(base, ".lsf.fits") ? g_strndup(base, strlen(base) - strlen(".lsf.fits")) : g_strdup(base);
	g_free(base);
	return stem;
}

/* With lsfm_newfile, Repair, Adopt and Compression write a new file, whose
 * name is asked here; otherwise they overwrite the file (*out NULL).
 * Returns FALSE when the user cancels. */
static gboolean ask_new_file(gchar **out) {
	*out = NULL;
	if (!info || !siril_toggle_get_active(w("lsfm_newfile")))
		return TRUE;
	gchar *stem = lsf_stem(info->file.path);
	gchar *gen = g_strrstr(stem, "_gen");	// a previous suggestion
	if (gen && gen[4] && strspn(gen + 4, "0123456789") == strlen(gen + 4))
		*gen = '\0';
	gchar *suggested = g_strdup_printf("%s_gen%d.lsf.fits", stem, info->file.gen + 1);
	gchar *file = choose_save_name(_("Living Stack files (*.lsf.fits)"), "*.lsf.fits", suggested);
	g_free(stem);
	g_free(suggested);
	if (!file)
		return FALSE;
	if (!g_str_has_suffix(file, ".lsf.fits")) {
		gchar *named = g_strdup_printf("%s.lsf.fits", file);
		g_free(file);
		file = named;
	}
	*out = file;
	return TRUE;
}

/* an lsf_util job that rewrites the file, or writes a new one */
static struct lsf_cmd_args *new_rewrite_args() {
	gchar *out;
	if (!info || processing_is_job_active() || !ask_new_file(&out)) {
		if (info && processing_is_job_active())
			PRINT_ANOTHER_THREAD_RUNNING;
		return NULL;
	}
	struct lsf_cmd_args *a = new_util_args();
	if (!a) {
		g_free(out);
		return NULL;
	}
	if (out && lsf_same_file(out, a->path))
		g_clear_pointer(&out, g_free);	// the file itself
	a->out = out;
	a->overwrite = out != NULL;	// the file chooser asked before replacing a file
	show_new_file = out != NULL;
	return a;
}

void on_lsfm_check_clicked(GtkButton *button, gpointer user_data) {
	struct lsf_cmd_args *a = new_util_args();
	if (!a)
		return;
	a->check = TRUE;
	start_job(a);
}

void on_lsfm_repair_clicked(GtkButton *button, gpointer user_data) {
	if (!siril_confirm_dialog(_("Repair the stack?"), _("The stacked image of the file is rebuilt from the "
					"accumulator, which undoes the edit made outside Siril."), _("Repair")))
		return;
	struct lsf_cmd_args *a = new_rewrite_args();
	if (!a)
		return;
	a->repair = TRUE;
	start_job(a);
}

void on_lsfm_adopt_clicked(GtkButton *button, gpointer user_data) {
	if (!siril_confirm_dialog(_("Adopt the edit?"), _("The edited pixels become the new mean of the accumulator. "
					"Their variance is not changed, so later rejection treats them as if they had been "
					"stacked this way. This is meant for small cosmetic fixes; process an exported image "
					"for anything else."), _("Adopt")))
		return;
	struct lsf_cmd_args *a = new_rewrite_args();
	if (!a)
		return;
	a->adopt = TRUE;
	start_job(a);
}

void on_lsfm_compress_clicked(GtkButton *button, gpointer user_data) {
	struct lsf_cmd_args *a = new_rewrite_args();
	if (!a)
		return;
	/* lsfm_compression: none, gzip2, gzip1 */
	guint sel = gtk_drop_down_get_selected(GTK_DROP_DOWN(w("lsfm_compression")));
	a->compression = sel == 1 ? LSF_COMP_GZIP2 : sel == 2 ? LSF_COMP_GZIP1 : LSF_COMP_NONE;
	a->compression_set = TRUE;
	start_job(a);
}

/* a file name chosen to save, or NULL */
static gchar *choose_save_name(const char *title_filter, const char *pattern, const char *suggested) {
	SirilFileChooser *fc = siril_fc_save(GTK_WINDOW(w("lsf_manager_dialog")), GTK_FILE_CHOOSER_ACTION_SAVE);
	siril_fc_set_current_folder_path(fc, com.wd);
	siril_fc_add_filter_pattern(fc, title_filter, pattern, TRUE);
	siril_fc_set_current_name(fc, suggested);
	gchar *file = siril_fc_run(fc) == GTK_RESPONSE_ACCEPT ? siril_fc_get_filename(fc) : NULL;
	siril_fc_destroy(fc);
	return file;
}

void on_lsfm_copy_clicked(GtkButton *button, gpointer user_data) {
	if (!info)
		return;
	gchar *stem = lsf_stem(info->file.path);
	gchar *suggested = g_strdup_printf("%s_copy.lsf.fits", stem);
	gchar *out = choose_save_name(_("Living Stack files (*.lsf.fits)"), "*.lsf.fits", suggested);
	g_free(stem);
	g_free(suggested);
	if (!out)
		return;
	if (!g_str_has_suffix(out, ".lsf.fits")) {
		gchar *named = g_strdup_printf("%s.lsf.fits", out);
		g_free(out);
		out = named;
	}
	struct lsf_cmd_args *a = new_util_args();
	if (!a) {
		g_free(out);
		return;
	}
	a->out = out;
	a->overwrite = TRUE;	// the file chooser asked before replacing a file
	start_job(a);
}

void on_lsfm_export_clicked(GtkButton *button, gpointer user_data) {
	if (!info || processing_is_job_active())
		return;
	gchar *base = g_path_get_basename(info->file.path);
	gchar *stem = g_str_has_suffix(base, ".lsf.fits") ? g_strndup(base, strlen(base) - strlen(".lsf.fits")) : g_strdup(base);
	gchar *suggested = g_strdup_printf("%s%s", stem, com.pref.ext);
	gchar *out = choose_save_name(_("FITS files"), "*.fit;*.fits;*.fts;*.fit.fz;*.fits.fz", suggested);
	g_free(base);
	g_free(stem);
	g_free(suggested);
	if (!out)
		return;
	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_EXPORT);
	if (!a) {
		g_free(out);
		return;
	}
	a->path = g_strdup(info->file.path);
	a->out = out;
	a->end_func = end_export;
	g_free(pending_out);
	pending_out = set_right_extension(out);
	gtk_label_set_text(GTK_LABEL(w("lsfm_result")), _("Working…"));
	start_job(a);
}

/* lsf_util -wcs= with a plate-solved image of the stack or a .wcs file */
void on_lsfm_import_wcs_clicked(GtkButton *button, gpointer user_data) {
	if (!info || processing_is_job_active())
		return;
	SirilFileChooser *fc = siril_fc_open(GTK_WINDOW(w("lsf_manager_dialog")), GTK_FILE_CHOOSER_ACTION_OPEN);
	siril_fc_set_current_folder_path(fc, com.wd);
	siril_fc_set_select_multiple(fc, FALSE);
	siril_fc_add_filter_pattern(fc, _("Plate-solved images and WCS files"),
			"*.fit;*.fits;*.fts;*.fit.fz;*.fits.fz;*.wcs;*.FIT;*.FITS;*.FTS;*.WCS", TRUE);
	gchar *file = siril_fc_run(fc) == GTK_RESPONSE_ACCEPT ? siril_fc_get_filename(fc) : NULL;
	siril_fc_destroy(fc);
	if (!file)
		return;
	struct lsf_cmd_args *a = new_rewrite_args();
	if (!a) {
		g_free(file);
		return;
	}
	a->wcs_file = file;
	start_job(a);
}

/* ── a WCS solved on a view of the file ────────────────────────────────── */

static gchar *view_wcs_file = NULL;	// the .wcs file written for lsf_util -wcs=

/* Siril removes the solutions it reads back (solve-field's .wcs) unless the
 * preferences keep them: this .wcs file is such an intermediate, since the
 * WCS is now in the Living Stack. It is kept when the import failed. */
static void end_view_wcs(const lsf_job_result *res) {
	set_cursor_waiting(FALSE);
	lsf_gui_report_write(res);
	if (!res->retval && view_wcs_file && !com.pref.astrometry.keep_wcs_files)
		g_unlink(view_wcs_file);
	g_clear_pointer(&view_wcs_file, g_free);
	if (res->retval) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("The WCS was not written into the Living Stack"),
				res->message ? res->message : _("See the log."));
		return;
	}
	/* the view shows the file as it is now, with the WCS it stores */
	open_single_image(res->path);
	siril_message_dialog(GTK_MESSAGE_INFO, _("WCS written into the Living Stack"),
			res->message ? res->message : _("Done."));
}

static gboolean read_view_info(lsf_info *in) {
	if (!single_image_is_loaded() || !com.uniq || !com.uniq->is_lsf || !com.uniq->filename)
		return FALSE;
	if (lsf_read_info(com.uniq->filename, in)) {
		lsf_info_free(in);
		return FALSE;
	}
	return TRUE;
}

/* the image shown is a view of an undistorted Living Stack, whose WCS must
 * stay linear: the astrometry dialog then solves at order 1 */
gboolean lsf_gui_view_is_undistorted(void) {
	lsf_info in;
	if (!read_view_info(&in))
		return FALSE;
	gboolean undistorted = in.file.settings.undistorted;
	lsf_info_free(&in);
	return undistorted;
}

/* After plate solving a view of a Living Stack (the astrometry dialog):
 * offers to write the solution into the file, through lsf_util -wcs= and a
 * .wcs file in the working directory. The solution is valid for the file
 * only if the view kept its geometry: the same LSFUID, size and row order. */
void lsf_gui_offer_view_wcs(void) {
	lsf_info in;
	if (!single_image_is_loaded() || !com.uniq || !com.uniq->is_lsf || !has_wcs(gfit))
		return;
	if (!read_view_info(&in)) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("The WCS cannot be written into the Living Stack"),
				_("The file of this view cannot be read as a Living Stack any more."));
		return;
	}
	const char *why = NULL;
	if (strcmp(in.file.uid, com.uniq->lsf_uid))
		why = _("The file was replaced by another Living Stack since this view was opened.");
	else if ((int) gfit->rx != in.file.rx || (int) gfit->ry != in.file.ry)
		why = _("The view and the Living Stack no longer have the same size: the view was cropped or resized, "
				"or frames grew the grid since it was opened. Open the file again and solve it.");
	else if (strcmp(gfit->keywords.row_order, in.file.primary.keywords.row_order))
		why = _("The view was flipped since it was opened. Open the file again and solve it.");
	if (why) {
		siril_message_dialog(GTK_MESSAGE_WARNING, _("The WCS cannot be written into the Living Stack"), (char *) why);
		lsf_info_free(&in);
		return;
	}
	gchar *base = g_path_get_basename(in.file.path);
	gchar *msg = g_strdup_printf(_("Write this solution into %s? The stacked image does not change; the file gets a new "
				"generation, and the view is opened again from it. The solution must be of the view as it was "
				"opened: not rotated, mirrored or otherwise transformed.%s"), base,
			in.file.settings.undistorted ? _(" The grid is undistorted: distortion terms are not kept.") : "");
	gboolean write = siril_confirm_dialog(_("Write the WCS into the Living Stack"), msg, _("Write"));
	g_free(msg);
	g_free(base);
	if (!write || processing_is_job_active()) {
		lsf_info_free(&in);
		return;
	}
	gchar *stem = lsf_stem(in.file.path), *name = g_strdup_printf("%s.lsf.wcs", stem);
	g_free(view_wcs_file);
	view_wcs_file = g_build_filename(com.wd, name, NULL);
	g_free(stem);
	g_free(name);
	struct lsf_cmd_args *a = NULL;
	if (save_wcs_fits(gfit, view_wcs_file) || !(a = lsf_cmd_args_new(LSF_CMD_UTIL))) {
		siril_message_dialog(GTK_MESSAGE_ERROR, _("The WCS was not written into the Living Stack"),
				_("The solution could not be saved to a .wcs file in the working directory."));
		g_clear_pointer(&view_wcs_file, g_free);
		lsf_info_free(&in);
		return;
	}
	a->path = g_strdup(in.file.path);
	a->wcs_file = g_strdup(view_wcs_file);
	a->end_func = end_view_wcs;
	lsf_info_free(&in);
	start_job(a);
}

void on_lsfm_maps_save_clicked(GtkButton *button, gpointer user_data) {
	if (!info || processing_is_job_active())
		return;
	gchar *prefix = g_strstrip(g_strdup(gtk_editable_get_text(GTK_EDITABLE(w("lsfm_maps_prefix")))));
	guint maps = 0;
	for (int i = 0; i < (int) G_N_ELEMENTS(map_kinds); i++) {
		GtkWidget *check = w(map_kinds[i].id);
		if (gtk_widget_get_sensitive(check) && siril_toggle_get_active(check))
			maps |= map_kinds[i].flag;
	}
	if (!prefix[0] || !maps) {
		gtk_label_set_text(GTK_LABEL(w("lsfm_maps_result")), _("Choose at least one map and a prefix."));
		g_free(prefix);
		return;
	}
	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_MAPS);
	if (!a) {
		g_free(prefix);
		return;
	}
	a->path = g_strdup(info->file.path);
	a->out = g_strdup(prefix);
	a->maps = maps;
	a->end_func = end_maps;
	g_free(pending_out);
	pending_out = prefix;
	gtk_label_set_text(GTK_LABEL(w("lsfm_maps_result")), _("Working…"));
	start_job(a);
}

/* ── window ────────────────────────────────────────────────────────────── */

void on_lsfm_close_clicked(GtkButton *button, gpointer user_data) {
	siril_close_dialog("lsf_manager_dialog");
}

gboolean lsf_manager_hide_on_delete(GtkWidget *widget) {
	siril_close_dialog("lsf_manager_dialog");
	return TRUE;
}

/* opens the dialog on path, or on the Living Stack loaded as the current
 * image, or else on the file shown before */
void lsf_manager_open(const gchar *path) {
	build_ledger_view();
	GtkEditable *entry = GTK_EDITABLE(w("lsfm_entry"));
	if (!path && single_image_is_loaded() && com.uniq->is_lsf)
		path = com.uniq->filename;
	if (path && g_strcmp0(path, gtk_editable_get_text(entry)))
		gtk_editable_set_text(entry, path);	// loaded by the changed handler
	else load_file();
	siril_open_dialog("lsf_manager_dialog");
}
