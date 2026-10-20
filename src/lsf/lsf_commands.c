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

/* LSF commands: lsf_init, lsf_append, lsf_export (FITS profile v0.4,
 * Sections 1 and 10). These are the only paths that write an LSF file. */

#include <string.h>
#include <glib/gstdio.h>

#include "yyjson.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "algos/siril_wcs.h"
#include "registration/registration.h"
#include "lsf/lsf.h"

typedef enum { LSF_CMD_INIT, LSF_CMD_APPEND, LSF_CMD_EXPORT } lsf_cmd_mode;

struct lsf_cmd_args {
	lsf_cmd_mode mode;
	gchar *seqname;
	gchar *path;
	gchar *out;		// export
	lsf_settings settings;	// init: requested settings
	gboolean hist_low_set, hist_high_set, hist_margin_set;
	double hist_low, hist_high, hist_margin;
	gboolean noreg, backup, overwrite;
};

static void free_cmd_args(struct lsf_cmd_args *a) {
	g_free(a->seqname);
	g_free(a->path);
	g_free(a->out);
	lsf_settings_free(&a->settings);
	free(a);
}

static gchar *findstar_to_json(const star_finder_params *sf) {
	gchar buf[6][G_ASCII_DTOSTR_BUF_SIZE];
	return g_strdup_printf("{\"radius\":%d,\"sigma\":%s,\"roundness\":%s,\"convergence\":%d,"
			"\"relax_checks\":%s,\"profile\":%d,\"min_beta\":%s,\"min_A\":%s,\"max_A\":%s,\"max_r\":%s}",
			sf->radius, g_ascii_formatd(buf[0], sizeof buf[0], "%.17g", sf->sigma),
			g_ascii_formatd(buf[1], sizeof buf[1], "%.17g", sf->roundness), sf->convergence,
			sf->relax_checks ? "true" : "false", (int) sf->profile,
			g_ascii_formatd(buf[2], sizeof buf[2], "%.17g", sf->min_beta),
			g_ascii_formatd(buf[3], sizeof buf[3], "%.17g", sf->min_A),
			g_ascii_formatd(buf[4], sizeof buf[4], "%.17g", sf->max_A),
			g_ascii_formatd(buf[5], sizeof buf[5], "%.17g", sf->max_r));
}

/* applies the star detection parameters stored in align.findstar */
static void findstar_from_json(const char *json, star_finder_params *sf) {
	if (!json) return;
	yyjson_doc *doc = yyjson_read(json, strlen(json), 0);
	if (!doc) return;
	yyjson_val *root = yyjson_doc_get_root(doc), *v;
	if ((v = yyjson_obj_get(root, "radius")) && yyjson_is_num(v)) sf->radius = (int) yyjson_get_num(v);
	if ((v = yyjson_obj_get(root, "sigma")) && yyjson_is_num(v)) sf->sigma = yyjson_get_num(v);
	if ((v = yyjson_obj_get(root, "roundness")) && yyjson_is_num(v)) sf->roundness = yyjson_get_num(v);
	if ((v = yyjson_obj_get(root, "convergence")) && yyjson_is_num(v)) sf->convergence = (int) yyjson_get_num(v);
	if ((v = yyjson_obj_get(root, "relax_checks")) && yyjson_is_bool(v)) sf->relax_checks = yyjson_get_bool(v);
	if ((v = yyjson_obj_get(root, "profile")) && yyjson_is_num(v)) sf->profile = (starprofile) (int) yyjson_get_num(v);
	if ((v = yyjson_obj_get(root, "min_beta")) && yyjson_is_num(v)) sf->min_beta = yyjson_get_num(v);
	if ((v = yyjson_obj_get(root, "min_A")) && yyjson_is_num(v)) sf->min_A = yyjson_get_num(v);
	if ((v = yyjson_obj_get(root, "max_A")) && yyjson_is_num(v)) sf->max_A = yyjson_get_num(v);
	if ((v = yyjson_obj_get(root, "max_r")) && yyjson_is_num(v)) sf->max_r = yyjson_get_num(v);
	yyjson_doc_free(doc);
}

static gchar *now_utc_iso(void) {
	GDateTime *now = g_date_time_new_now_utc();
	gchar *s = g_date_time_format(now, "%Y-%m-%dT%H:%M:%S");
	g_date_time_unref(now);
	return s;
}

static sequence *load_batch(const char *name) {
	sequence *seq = load_sequence(name, NULL);
	if (!seq)
		return NULL;
	if (seq_check_basic_data(seq, FALSE) < 0) {
		free_sequence(seq, TRUE);
		return NULL;
	}
	return seq;
}

static void apply_hist_overrides(struct lsf_cmd_args *a, lsf_settings *s) {
	if (a->hist_low_set) s->hist_low = a->hist_low;
	if (a->hist_high_set) s->hist_high = a->hist_high;
	if (a->hist_margin_set) s->hist_margin = a->hist_margin;
}

static int lsf_init_run(struct lsf_cmd_args *a) {
	lsf_settings *s = &a->settings;
	lsf_history h = { 0 };
	fits primary = { 0 };
	GPtrArray *rows = g_ptr_array_new_with_free_func(lsf_ledger_row_free);
	int retval = 1;
	gchar *uid = NULL;

	if (g_file_test(a->path, G_FILE_TEST_EXISTS) && !a->overwrite) {
		siril_log_error(_("LSF: %s already exists, use lsf_append to add frames to it, or -overwrite\n"), a->path);
		g_ptr_array_free(rows, TRUE);
		return 1;
	}
	sequence *seq = load_batch(a->seqname);
	if (!seq) {
		g_ptr_array_free(rows, TRUE);
		return 1;
	}
	int reglayer = get_registration_layer(seq);
	if (layer_has_usable_registration(seq, reglayer)) {
		siril_log_error(_("LSF: the sequence %s has registration data that is not applied. lsf_init needs "
					"frames already on a common grid: apply the registration first (seqapplyreg or "
					"register without -2pass)\n"), a->seqname);
		goto end;
	}

	/* sensor and fixed settings */
	fits refmeta = { 0 };
	int ref = sequence_find_refimage(seq);
	if (!seq_read_frame_metadata(seq, ref, &refmeta)) {
		s->is_cfa = seq->nb_layers == 3;
		const char *pat = refmeta.keywords.bayer_pattern;
		if (s->is_cfa)
			g_strlcpy(s->bayer, (!strcmp(pat, "RGGB") || !strcmp(pat, "GRBG") || !strcmp(pat, "GBRG")
						|| !strcmp(pat, "BGGR")) ? pat : "UNKNOWN", sizeof s->bayer);
		s->binning = refmeta.keywords.binning_x ? (int) refmeta.keywords.binning_x : 1;
		g_free(s->filter);
		s->filter = g_strdup(refmeta.keywords.filter);
		clearfits(&refmeta);
	}
	if (s->rgb_equal && seq->nb_layers != 3)
		s->rgb_equal = FALSE;
	g_free(s->findstar_json);
	s->findstar_json = findstar_to_json(&com.pref.starfinder_conf);
	g_free(s->created_utc);
	s->created_utc = now_utc_iso();
	if (!a->hist_low_set)
		s->hist_low = lsf_rejection_has_sigma(s->rej) ? s->sig[0] : 3.0;
	if (!a->hist_high_set)
		s->hist_high = lsf_rejection_has_sigma(s->rej) ? s->sig[1] : 3.0;
	apply_hist_overrides(a, s);

	struct lsf_batch_args bargs = { 0 };
	bargs.seq = seq;
	bargs.settings = s;
	bargs.hist = &h;
	bargs.init = TRUE;
	bargs.ref_image = ref;
	bargs.batch = 1;
	bargs.gen = 1;
	bargs.prev_frames = 0;
	bargs.seq_name = seq->seqname;
	bargs.rows = rows;
	bargs.primary = &primary;
	if (lsf_process_batch(&bargs))
		goto end;
	if (lsf_build_primary(&primary, s, &h))
		goto end;
	uid = g_uuid_string_random();
	if (lsf_write_file(a->path, NULL, s, &h, &primary, uid, 1, rows, 1, FALSE))
		goto end;
	siril_log_info(_("LSF file %s created: %u frames, generation 1\n"),
			a->path, rows->len);
	retval = 0;
end:
	g_free(uid);
	clearfits(&primary);
	lsf_history_free(&h);
	g_ptr_array_free(rows, TRUE);
	free_sequence(seq, TRUE);
	return retval;
}

/* registers the batch against HDU 0 of the LSF file with Siril's global star
 * alignment and an external reference (Section 10.5). Returns the registered
 * sequence, written with the lsf_r_ prefix. */
#define LSF_REG_PREFIX "lsf_r_"
static sequence *register_against_stack(sequence *seq, const char *lsf_path, const lsf_settings *s) {
	star_finder_params backup = com.pref.starfinder_conf;
	findstar_from_json(s->findstar_json, &com.pref.starfinder_conf);

	struct registration_args *regargs = calloc(1, sizeof(struct registration_args));
	regargs->func = register_star_alignment;
	regargs->seq = seq;
	regargs->reference_image = sequence_find_refimage(seq);
	regargs->layer = seq->nb_layers == 3 ? 1 : 0;
	regargs->filters.filter_included = TRUE;
	regargs->run_in_thread = FALSE;
	regargs->no_output = FALSE;
	regargs->output_scale = 1.f;
	regargs->prefix = strdup(LSF_REG_PREFIX);
	regargs->min_pairs = s->min_pairs;
	regargs->max_stars_candidates = MAX_STARS_FITTED;
	regargs->type = s->transform;
	regargs->interpolation = s->interpolation;
	regargs->clamp = s->interpolation == OPENCV_CUBIC || s->interpolation == OPENCV_LANCZOS4;
	regargs->undistort = DISTO_UNDEF;
	regargs->load_new_sequence = FALSE;
	regargs->external_ref_path = g_strdup(lsf_path);
	regargs->use_external_ref = TRUE;

	siril_log_message(_("LSF: registering %d frames against the Stack of %s\n"), seq->selnum, lsf_path);
	int ret = register_star_alignment(regargs);
	com.pref.starfinder_conf = backup;
	if (!ret)
		writeseqfile(seq);

	free(regargs->prefix);
	g_free(regargs->external_ref_path);
	g_free(regargs->new_seq_name);
	if (regargs->wcsref)
		wcsfree(regargs->wcsref);
	free(regargs->sfargs);
	free(regargs);
	if (ret) {
		siril_log_error(_("LSF: registration against the Stack failed\n"));
		return NULL;
	}
	gchar *rname = g_strdup_printf("%s%s", LSF_REG_PREFIX, seq->seqname);
	sequence *rseq = load_batch(rname);
	g_free(rname);
	return rseq;
}

static int lsf_append_run(struct lsf_cmd_args *a) {
	lsf_file f;
	int retval = 1;
	sequence *seq = NULL, *rseq = NULL;
	int *src_index = NULL;
	GPtrArray *rows = g_ptr_array_new_with_free_func(lsf_ledger_row_free);

	if (lsf_open_and_validate(a->path, &f)) {
		lsf_file_free(&f);
		g_ptr_array_free(rows, TRUE);
		return 1;
	}
	apply_hist_overrides(a, &f.settings);
	seq = load_batch(a->seqname);
	if (!seq)
		goto end;
	if (seq->nb_layers != f.nchan) {
		siril_log_error(_("LSF: the file has %d channel(s), the batch has %d\n"), f.nchan, seq->nb_layers);
		goto end;
	}

	sequence *batch = seq;
	if (!a->noreg) {
		rseq = register_against_stack(seq, a->path, &f.settings);
		if (!rseq)
			goto end;
		/* map registered frames back to their sources */
		src_index = malloc(rseq->number * sizeof(int));
		for (int i = 0; i < rseq->number; i++) {
			src_index[i] = -1;
			for (int j = 0; j < seq->number; j++)
				if (seq->imgparam[j].filenum == rseq->imgparam[i].filenum) {
					src_index[i] = j;
					break;
				}
			if (src_index[i] < 0) {
				siril_log_error(_("LSF: cannot map registered frame %d back to its source\n"), rseq->imgparam[i].filenum);
				goto end;
			}
		}
		int dropped = seq->selnum - rseq->number;
		if (dropped > 0)
			siril_log_warning(_("LSF: %d frame(s) could not be registered against the Stack and are not merged\n"), dropped);
		batch = rseq;
	} else {
		int reglayer = get_registration_layer(seq);
		if (layer_has_usable_registration(seq, reglayer)) {
			siril_log_error(_("LSF: -noreg needs frames already on the LSF grid, but %s has registration data that is not applied\n"),
					a->seqname);
			goto end;
		}
	}

	int gen = f.gen + 1, batchnum = f.nb_batches + 1;
	struct lsf_batch_args bargs = { 0 };
	bargs.seq = batch;
	bargs.src_seq = a->noreg ? NULL : seq;
	bargs.src_index = src_index;
	bargs.settings = &f.settings;
	bargs.hist = &f.hist;
	bargs.init = FALSE;
	bargs.batch = batchnum;
	bargs.gen = gen;
	bargs.prev_frames = f.ledger_rows;
	bargs.seq_name = seq->seqname;
	bargs.rows = rows;
	bargs.primary = &f.primary;
	if (lsf_process_batch(&bargs))
		goto end;
	if (lsf_build_primary(&f.primary, &f.settings, &f.hist))
		goto end;
	if (lsf_write_file(a->path, &f, &f.settings, &f.hist, &f.primary, f.uid, gen, rows, batchnum, a->backup))
		goto end;
	siril_log_info(_("LSF: %u frames appended to %s, %ld frames in total, generation %d\n"),
			rows->len, a->path, f.ledger_rows + rows->len, gen);
	retval = 0;
end:
	if (rseq) {
		remove_prefixed_sequence_files(seq, LSF_REG_PREFIX);
		free_sequence(rseq, TRUE);
	}
	free(src_index);
	if (seq)
		free_sequence(seq, TRUE);
	g_ptr_array_free(rows, TRUE);
	lsf_file_free(&f);
	return retval;
}

static gpointer lsf_worker(gpointer p) {
	struct lsf_cmd_args *a = p;
	int retval;
	switch (a->mode) {
	case LSF_CMD_INIT:
		retval = lsf_init_run(a);
		break;
	case LSF_CMD_APPEND:
		retval = lsf_append_run(a);
		break;
	default:
		retval = lsf_export_master(a->path, a->out);
		if (!retval)
			siril_log_message(_("LSF: master of %s exported to %s\n"), a->path, a->out);
		break;
	}
	free_cmd_args(a);
	siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

static gboolean parse_double_opt(const char *word_i, const char *prefix, double *out, gboolean *set) {
	if (!g_str_has_prefix(word_i, prefix))
		return FALSE;
	const char *value = word_i + strlen(prefix);
	gchar *end;
	double v = g_ascii_strtod(value, &end);
	if (end == value || v < 0.0) {
		siril_log_error(_("Wrong parameter values. %s is not valid, aborting.\n"), word_i);
		*set = -1;
		return TRUE;
	}
	*out = v;
	*set = TRUE;
	return TRUE;
}

/* options common to lsf_init and lsf_append; returns 1 if consumed, -1 on error */
static int parse_common_option(struct lsf_cmd_args *a, const char *w) {
	if (parse_double_opt(w, "-hist_low=", &a->hist_low, &a->hist_low_set))
		return a->hist_low_set == TRUE ? 1 : -1;
	if (parse_double_opt(w, "-hist_high=", &a->hist_high, &a->hist_high_set))
		return a->hist_high_set == TRUE ? 1 : -1;
	if (parse_double_opt(w, "-hist_margin=", &a->hist_margin, &a->hist_margin_set)) {
		if (a->hist_margin_set == TRUE && a->hist_margin < 1.0) {
			siril_log_error(_("-hist_margin must be at least 1\n"));
			return -1;
		}
		return a->hist_margin_set == TRUE ? 1 : -1;
	}
	return 0;
}

/* lsf_init sequencename file.lsf.fits [rej [type] sigma_low sigma_high] [-nonorm, -norm=] [-fastnorm]
 *   [-rgb_equal] [-weight=] [-output_norm] [-hist_low=] [-hist_high=] [-hist_margin=] [-overwrite] */
int process_lsf_init(int nb) {
	struct lsf_cmd_args *a = calloc(1, sizeof(struct lsf_cmd_args));
	a->mode = LSF_CMD_INIT;
	a->seqname = g_strdup(word[1]);
	a->path = g_strdup(word[2]);
	lsf_settings *s = &a->settings;
	lsf_settings_init(s);
	s->norm = NO_NORM;	// as the stack command
	s->rej = NO_REJEC;

	int i = 3;
	if (word[i] && (!strcmp(word[i], "rej") || !strcmp(word[i], "mean"))) {
		i++;
		int shift = 1;
		const char *t = word[i] ? word[i] : "";
		if (!strcmp(t, "p") || !strcmp(t, "percentile")) s->rej = PERCENTILE;
		else if (!strcmp(t, "s") || !strcmp(t, "sigma")) s->rej = SIGMA;
		else if (!strcmp(t, "a") || !strcmp(t, "mad")) s->rej = MAD;
		else if (!strcmp(t, "m") || !strcmp(t, "median")) s->rej = SIGMEDIAN;
		else if (!strcmp(t, "l") || !strcmp(t, "linear")) s->rej = LINEARFIT;
		else if (!strcmp(t, "w") || !strcmp(t, "winsorized")) s->rej = WINSORIZED;
		else if (!strcmp(t, "g") || !strcmp(t, "generalized")) s->rej = GESDT;
		else if (!strcmp(t, "n") || !strcmp(t, "none")) s->rej = NO_REJEC;
		else { s->rej = WINSORIZED; shift = 0; }
		i += shift;
		if (s->rej != NO_REJEC || (word[i] && string_is_a_number(word[i]))) {
			gchar *end1 = NULL, *end2 = NULL;
			if (!word[i] || !word[i + 1] || !string_is_a_number(word[i]) || !string_is_a_number(word[i + 1])) {
				siril_log_error(_("The average stacking with rejection requires two extra arguments: sigma low and high.\n"));
				goto failure;
			}
			s->sig[0] = g_ascii_strtod(word[i], &end1);
			s->sig[1] = g_ascii_strtod(word[i + 1], &end2);
			if (s->sig[0] < 0.0 || s->sig[1] < 0.0) {
				siril_log_error(_("The average stacking with rejection requires two extra arguments: sigma low and high.\n"));
				goto failure;
			}
			if ((s->rej == GESDT || s->rej == PERCENTILE) && (s->sig[0] > 1.0 || s->sig[1] > 1.0)) {
				siril_log_error(_("Extra parameters of this rejection algorithm must be between 0 and 1.\n"));
				goto failure;
			}
			i += 2;
		}
	}
	for (; i < nb && word[i]; i++) {
		const char *w = word[i];
		int c = parse_common_option(a, w);
		if (c < 0) goto failure;
		if (c > 0) continue;
		if (!strcmp(w, "-nonorm") || !strcmp(w, "-no_norm")) s->norm = NO_NORM;
		else if (g_str_has_prefix(w, "-norm=")) {
			if (!lsf_norm_from_name(w + 6, &s->norm)) {
				siril_log_error(_("Unknown argument to %s, aborting.\n"), w);
				goto failure;
			}
		}
		else if (!strcmp(w, "-fastnorm")) s->lite_norm = TRUE;
		else if (!strcmp(w, "-rgb_equal")) s->rgb_equal = TRUE;
		else if (!strcmp(w, "-output_norm")) s->outnorm = TRUE;
		else if (!strcmp(w, "-overwrite")) a->overwrite = TRUE;
		else if (g_str_has_prefix(w, "-weight=")) {
			const char *v = w + 8;
			if (!strcmp(v, "noise")) s->weighting = NOISE_WEIGHT;
			else if (!strcmp(v, "nbstars")) s->weighting = NBSTARS_WEIGHT;
			else if (!strcmp(v, "wfwhm")) s->weighting = WFWHM_WEIGHT;
			else if (!strcmp(v, "none")) s->weighting = NO_WEIGHT;
			else {
				siril_log_error(_("LSF: weighting %s is not supported (noise, wfwhm, nbstars)\n"), v);
				goto failure;
			}
		}
		else {
			siril_log_error(_("Unknown parameter %s, aborting.\n"), w);
			goto failure;
		}
	}
	if (s->norm == NO_NORM) {
		s->lite_norm = FALSE;
		s->rgb_equal = FALSE;
	}
	if (!start_in_new_thread(lsf_worker, a)) {
		free_cmd_args(a);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
failure:
	free_cmd_args(a);
	return CMD_ARG_ERROR;
}

/* lsf_append file.lsf.fits sequencename [-noreg] [-backup] [-hist_low=] [-hist_high=] [-hist_margin=] */
int process_lsf_append(int nb) {
	struct lsf_cmd_args *a = calloc(1, sizeof(struct lsf_cmd_args));
	a->mode = LSF_CMD_APPEND;
	a->path = g_strdup(word[1]);
	a->seqname = g_strdup(word[2]);
	for (int i = 3; i < nb && word[i]; i++) {
		const char *w = word[i];
		int c = parse_common_option(a, w);
		if (c < 0) goto failure;
		if (c > 0) continue;
		if (!strcmp(w, "-noreg")) a->noreg = TRUE;
		else if (!strcmp(w, "-backup")) a->backup = TRUE;
		else {
			siril_log_error(_("Unknown parameter %s, aborting.\n"), w);
			goto failure;
		}
	}
	if (!start_in_new_thread(lsf_worker, a)) {
		free_cmd_args(a);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
failure:
	free_cmd_args(a);
	return CMD_ARG_ERROR;
}

/* lsf_export file.lsf.fits output.fit */
int process_lsf_export(int nb) {
	struct lsf_cmd_args *a = calloc(1, sizeof(struct lsf_cmd_args));
	a->mode = LSF_CMD_EXPORT;
	a->path = g_strdup(word[1]);
	a->out = g_strdup(word[2]);
	if (!start_in_new_thread(lsf_worker, a)) {
		free_cmd_args(a);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}
