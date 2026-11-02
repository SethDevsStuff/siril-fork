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

/* LSF commands: lsf_init, lsf_append, lsf_export, lsf_util, lsf_maps
 * (FITS profile v0.4, Sections 1, 8 and 10). These are the only paths that
 * write an LSF file. */

#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>

#include <glib/gstdio.h>

#include "algos/astrometry_solver.h"
#include "algos/siril_wcs.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "core/proto.h"
#include "core/sequence_filtering.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "io/single_image.h"
#include "lsf/lsf.h"
#include "lsf/lsf_commands.h"
#include "lsf/lsf_fits.h"
#include "lsf/lsf_math.h"
#include "lsf/lsf_settings.h"
#include "lsf/lsf_stack.h"
#include "lsf/lsf_wcs.h"
#include "lsf/lsf_weights.h"
#include "opencv/opencv.h"
#include "registration/registration.h"

#include "yyjson.h"

/* the defaults of the commands, which the parsers then change */
struct lsf_cmd_args *lsf_cmd_args_new(lsf_cmd_mode mode) {
	struct lsf_cmd_args *a = calloc(1, sizeof(struct lsf_cmd_args));
	if (!a) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	a->mode = mode;
	a->framing = FRAMING_MAX;
	a->expand_step = LSF_EXPAND_STEP;
	a->max_threads = com.max_thread;
	lsf_settings_init(&a->settings);
	a->settings.norm = NO_NORM;	// as the stack command
	a->settings.rej = NO_REJEC;
	return a;
}

void lsf_cmd_args_free(struct lsf_cmd_args *a) {
	if (!a)
		return;
	g_free(a->seqname);
	g_free(a->path);
	g_free(a->out);
	g_free(a->export_out);
	g_free(a->message);
	g_free(a->image_indices);
	g_free(a->wcs_file);
	g_free(a->disto.filename);
	lsf_settings_free(&a->settings);
	free(a);
}

static void add_message(struct lsf_cmd_args *a, const char *fmt, ...);

/* keeps the result of lsf_write_file() for the GUI: a failed replace puts
 * what the user should do in the message of the job. Returns status. */
static int note_write(struct lsf_cmd_args *a, int status, gchar *replace_msg) {
	a->write_status = status;
	if (replace_msg) {
		add_message(a, "%s", replace_msg);
		g_free(replace_msg);
	}
	return status;
}

/* adds a line to the message of the job, which the GUI shows */
static void add_message(struct lsf_cmd_args *a, const char *fmt, ...) {
	va_list args;
	va_start(args, fmt);
	gchar *line = g_strdup_vprintf(fmt, args);
	va_end(args);
	gchar *msg = a->message ? g_strdup_printf("%s\n%s", a->message, line) : g_strdup(line);
	g_free(line);
	g_free(a->message);
	a->message = msg;
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

/* the rejection as the log and the HISTORY cards show it: the method, then
 * its two values unless it is none */
static gchar *rejection_text(rejection r, const double sig[2]) {
	if (r == NO_REJEC)
		return g_strdup(lsf_rejection_name(r));
	return g_strdup_printf("%s %g %g", lsf_rejection_name(r), sig[0], sig[1]);
}

/* The frames the GUI's stacking filters chose, in place of
 * the selection saved in the sequence file. The change stays in memory:
 * the caller keeps lsf_process_batch() from writing it to the .seq file
 * of the user. Returns TRUE when a list was applied. */
static gboolean apply_frame_list(sequence *seq, const struct lsf_cmd_args *a) {
	if (!a->image_indices || a->nb_images < 1)
		return FALSE;
	for (int i = 0; i < seq->number; i++)
		seq->imgparam[i].incl = FALSE;
	int n = 0;
	for (int k = 0; k < a->nb_images; k++) {
		int i = a->image_indices[k];
		if (i >= 0 && i < seq->number && !seq->imgparam[i].incl) {
			seq->imgparam[i].incl = TRUE;
			n++;
		}
	}
	seq->selnum = n;
	return TRUE;
}

static void apply_hist_overrides(struct lsf_cmd_args *a, lsf_settings *s) {
	if (a->hist_low_set) s->hist_low = a->hist_low;
	if (a->hist_high_set) s->hist_high = a->hist_high;
	if (a->hist_margin_set) s->hist_margin = a->hist_margin;
}

#define LSF_REG_PREFIX "lsf_r_"
static sequence *apply_init_registration(sequence *seq, const struct lsf_cmd_args *a, int *ref,
		int origin[2], Homography *base, gboolean *undistorted);
static int *map_to_sources(sequence *seq, sequence *rseq);

static gboolean any_pad(const int pad[4]) {
	return pad[LSF_PAD_LEFT] || pad[LSF_PAD_RIGHT] || pad[LSF_PAD_FIRST] || pad[LSF_PAD_LAST];
}

static const char *lsf_framing_name(framing_type f) {
	switch (f) {
	case FRAMING_CURRENT: return "current";
	case FRAMING_MIN: return "min";
	case FRAMING_COG: return "cog";
	default: return "max";
	}
}

static int lsf_init_run(struct lsf_cmd_args *a) {
	lsf_settings *s = &a->settings;
	lsf_history h = { 0 };
	lsf_wref wref = { 0 };
	fits primary = { 0 };
	GPtrArray *rows = g_ptr_array_new_with_free_func(lsf_ledger_row_free);
	int retval = 1;
	gchar *uid = NULL;
	sequence *rseq = NULL;
	int *src_index = NULL;
	int origin[2] = { 0, 0 };
	Homography base;

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
	gboolean listed = apply_frame_list(seq, a);
	/* registration computed but not applied: lsf_init applies it, by default
	 * onto the union of all the frames (Section 10.1) */
	int ref = sequence_find_refimage(seq);
	sequence *batch = seq;
	int batch_ref = ref;
	if (layer_has_usable_registration(seq, get_registration_layer(seq))) {
		/* Siril's apply undistorts the frames of a registration made with
		 * distortion correction: the grid is then undistorted, for good, and
		 * every append undistorts its frames with the same model (Section 10.4) */
		int rl = get_registration_layer(seq);
		disto_source ds = seq->distoparam && rl >= 0 ? seq->distoparam[rl].index : DISTO_UNDEF;
		if (ds == DISTO_IMAGE || ds == DISTO_FILE || ds == DISTO_MASTER) {
			if (!seq->distoparam[rl].filename && ds != DISTO_MASTER) {
				siril_log_error(_("LSF: the distortion model of the registration of %s is not named in its .seq file, "
							"register it again\n"), a->seqname);
				goto end;
			}
			s->undistorted = TRUE;
			lsf_disto_copy(&s->disto, &seq->distoparam[rl]);
			if (s->disto.index == DISTO_IMAGE)
				s->disto.index = DISTO_FILE;	// Siril saved the image's WCS to this file
			if (s->disto.index == DISTO_MASTER && !s->disto.filename)
				s->disto.filename = g_strdup(com.pref.prepro.disto_lib);
			/* the .seq file names the model relative to the working directory;
			 * the file keeps it whatever directory the appends run from */
			if (s->disto.index == DISTO_FILE) {
				gchar *abs = g_canonicalize_filename(s->disto.filename, NULL);
				g_free(s->disto.filename);
				s->disto.filename = abs;
			}
			gchar *model = lsf_disto_describe(&s->disto);
			siril_log_message(_("LSF: %s was registered with distortion correction (%s): the frames are undistorted, "
						"and so is the grid. Every append will undistort its frames with this model, or with "
						"another one given to lsf_append with -disto=\n"), a->seqname, model);
			g_free(model);
		}
		gboolean undistorted = FALSE;
		rseq = apply_init_registration(seq, a, &ref, origin, &base, &undistorted);
		if (!rseq || !(src_index = map_to_sources(seq, rseq)))
			goto end;
		/* an astrometric registration (seqplatesolve): Siril's apply
		 * undistorts each frame with the SIP terms of its own WCS, when it
		 * has some, and the appends do the same through their WCS */
		if (ds == DISTO_FILES) {
			if (undistorted) {
				s->undistorted = TRUE;
				g_free(s->disto.filename);
				s->disto = (disto_params) { .index = DISTO_FILES };
				siril_log_message(_("LSF: %s was registered through plate solving: each frame is undistorted with the "
							"SIP terms of its own WCS, and the grid is undistorted. Every append must be plate-solved "
							"and registered through the WCS (lsf_append -wcsreg), or given one model with -disto=\n"),
						a->seqname);
			} else siril_log_message(_("LSF: %s was registered through plate solving, without distortion terms: the "
						"frames are not undistorted\n"), a->seqname);
		}
		batch = rseq;
		batch_ref = 0;
		for (int i = 0; i < rseq->number; i++)
			if (src_index[i] == ref)
				batch_ref = i;
		s->origin_x = origin[0];
		s->origin_y = origin[1];
	} else if (a->framing_set) {
		siril_log_warning(_("LSF: the frames of %s are already on a common grid, -framing is ignored\n"), a->seqname);
	}

	/* sensor and fixed settings */
	fits refmeta = { 0 };
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
	s->rejmaps = a->rejmaps;
	g_free(s->findstar_json);
	s->findstar_json = findstar_to_json(&com.pref.starfinder_conf);
	g_free(s->created_utc);
	s->created_utc = lsf_now_utc_iso();
	if (!a->hist_low_set)
		s->hist_low = lsf_rejection_has_sigma(s->rej) && s->sig[0] > 0.0 ? s->sig[0] : 3.0;
	if (!a->hist_high_set)
		s->hist_high = lsf_rejection_has_sigma(s->rej) && s->sig[1] > 0.0 ? s->sig[1] : 3.0;
	apply_hist_overrides(a, s);

	struct lsf_batch_args bargs = { 0 };
	bargs.seq = batch;
	bargs.src_seq = rseq ? seq : NULL;
	bargs.src_index = src_index;
	bargs.xfrm_base = rseq ? &base : NULL;
	bargs.framing = lsf_framing_name(a->framing);
	bargs.disto = s->undistorted ? &s->disto : NULL;
	bargs.settings = s;
	bargs.hist = &h;
	bargs.init = TRUE;
	bargs.ref_image = batch_ref;
	bargs.batch = 1;
	bargs.gen = 1;
	bargs.prev_frames = 0;
	bargs.seq_name = seq->seqname;
	bargs.rows = rows;
	bargs.primary = &primary;
	bargs.wref = &wref;
	bargs.keep_seq_file = listed && batch == seq;
	bargs.max_threads = a->max_threads;
	bargs.feather = a->feather;
	bargs.overlap_norm = a->overlap_norm;
	if (lsf_process_batch(&bargs))
		goto end;
	/* a HISTORY card holds 72 characters: the undistortion has its own */
	if (rseq)
		primary.history = g_slist_append(primary.history, g_strdup_printf(
					"LSF gen 1: registration applied, %s framing, origin %d %d", bargs.framing,
					origin[0], origin[1]));
	if (rseq && s->undistorted)
		primary.history = g_slist_append(primary.history, g_strdup(s->disto.index == DISTO_FILES ?
					"LSF gen 1: grid undistorted, each frame with its own SIP terms" :
					"LSF gen 1: grid undistorted"));
	if (lsf_build_primary(&primary, s, &h))
		goto end;
	uid = g_uuid_string_random();
	gchar *replace_msg = NULL;
	int written = lsf_write_file(a->path, NULL, s, &h, &wref, &primary, uid, 1, rows, 1, FALSE, a->compression, &replace_msg);
	if (note_write(a, written, replace_msg))
		goto end;
	siril_log_info(_("LSF file %s created: %u frames, generation 1\n"),
			a->path, rows->len);
	retval = 0;
end:
	g_free(uid);
	clearfits(&primary);
	lsf_history_free(&h);
	lsf_wref_free(&wref);
	g_ptr_array_free(rows, TRUE);
	if (rseq) {
		remove_prefixed_sequence_files(seq, LSF_REG_PREFIX);
		free_sequence(rseq, TRUE);
	}
	free(src_index);
	free_sequence(seq, TRUE);
	return retval;
}

/* registration arguments for Siril's star alignment and apply steps, or NULL */
static struct registration_args *new_regargs(sequence *seq, const lsf_settings *s) {
	struct registration_args *regargs = calloc(1, sizeof(struct registration_args));
	if (!regargs) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	regargs->prefix = strdup(LSF_REG_PREFIX);
	if (!regargs->prefix) {
		PRINT_ALLOC_ERR;
		free(regargs);
		return NULL;
	}
	regargs->seq = seq;
	regargs->reference_image = sequence_find_refimage(seq);
	regargs->layer = seq->nb_layers == 3 ? 1 : 0;
	regargs->filters.filter_included = TRUE;
	regargs->run_in_thread = FALSE;
	regargs->output_scale = 1.f;
	regargs->type = s->transform;
	regargs->interpolation = s->interpolation;
	regargs->clamp = s->interpolation == OPENCV_CUBIC || s->interpolation == OPENCV_LANCZOS4;
	regargs->undistort = DISTO_UNDEF;
	regargs->framing = FRAMING_CURRENT;
	regargs->load_new_sequence = FALSE;
	return regargs;
}

static void free_regargs(struct registration_args *regargs) {
	free(regargs->prefix);
	g_free(regargs->external_ref_path);
	g_free(regargs->new_seq_name);
	if (regargs->wcsref)
		wcsfree(regargs->wcsref);
	free(regargs->sfargs);
	free(regargs);
}

/* grid expansion (Section 10.5): the padding of each side needed by the
 * footprint box of the registered frames, or 0 when no expansion is wanted
 * or needed. Fails when the expansion is refused. */
static int plan_expansion(const lsf_file *f, const struct lsf_cmd_args *a, const double box[4], int pad[4]) {
	int need[4];
	memset(pad, 0, 4 * sizeof(int));
	if (!lsf_grid_overhang(box, f->rx, f->ry, need))
		return 0;
	if (a->framing == FRAMING_CURRENT) {
		siril_log_warning(_("LSF: the batch reaches beyond the grid by %d px on the left, %d on the right, "
					"%d before row 1 and %d after the last row; these pixels are cropped (-framing=current)\n"),
				need[LSF_PAD_LEFT], need[LSF_PAD_RIGHT], need[LSF_PAD_FIRST], need[LSF_PAD_LAST]);
		return 0;
	}
	lsf_expand_pads(need, a->expand_step, pad);
	int nrx = f->rx + pad[LSF_PAD_LEFT] + pad[LSF_PAD_RIGHT];
	int nry = f->ry + pad[LSF_PAD_FIRST] + pad[LSF_PAD_LAST];
	if (!a->expand_force && (nrx > 2 * f->rx || nry > 2 * f->ry)) {
		siril_log_error(_("LSF: this batch would grow the grid from %dx%d to %dx%d, more than twice on one axis, "
					"which usually means a wrong registration. Check the registration, or use -expand_force "
					"to grow the grid anyway, as a mosaic panel appended beside the Stack may need, or "
					"-framing=current to crop the batch to the grid\n"),
				f->rx, f->ry, nrx, nry);
		return 1;
	}
	guint64 bytes = lsf_history_bytes((guint64) nrx * nry * f->nchan, f->settings.rejmaps, f->has_exp);
	int max_mb = get_max_memory_in_MB();
	if (bytes > (guint64) max_mb * BYTES_IN_A_MB) {
		siril_log_error(_("LSF: the accumulator of the expanded %dx%d grid needs %.0f MB, more than the %d MB "
					"Siril may use. Use -framing=current to crop the batch to the grid, or raise the memory limit\n"),
				nrx, nry, (double) bytes / BYTES_IN_A_MB, max_mb);
		return 1;
	}
	siril_log_message(_("LSF: the grid grows from %dx%d to %dx%d: %d px on the left, %d on the right, "
				"%d before row 1, %d after the last row\n"), f->rx, f->ry, nrx, nry,
			pad[LSF_PAD_LEFT], pad[LSF_PAD_RIGHT], pad[LSF_PAD_FIRST], pad[LSF_PAD_LAST]);
	return 0;
}

/* first pass of the registration against HDU 0 of the LSF file (Section
 * 10.4): Siril's global star alignment with the file as external
 * reference, with the star detection, transform and pairs stored in the
 * file. The transforms stay in memory, in seq->regparam[*layer]: the .seq
 * file of the batch and its cached star lists are left as the user made
 * them. A registration made beforehand replaces this pass. */
static int compute_against_stack(sequence *seq, const lsf_file *f, const disto_params *model, int *layer) {
	const lsf_settings *s = &f->settings;
	/* the star detection settings stored in the file, over a copy of the
	 * preferences for what they do not hold; the preferences are not changed */
	star_finder_params sf = com.pref.starfinder_conf;
	lsf_findstar_from_json(s->findstar_json, &sf);

	struct registration_args *regargs = new_regargs(seq, s);
	if (!regargs)
		return 1;
	regargs->func = register_star_alignment;
	regargs->no_output = TRUE;
	regargs->no_starlist = TRUE;
	regargs->sf_params = &sf;
	regargs->min_pairs = s->min_pairs;
	regargs->max_stars_candidates = MAX_STARS_FITTED;
	regargs->external_ref_path = g_strdup(f->path);
	regargs->use_external_ref = TRUE;
	/* an undistorted grid: the stars of the frames are undistorted with the
	 * model before matching, the Stack's are not, since it is undistorted
	 * already. Siril then records the model on the registration layer, and
	 * its apply undistorts the frames with it. */
	if (model) {
		regargs->undistort = model->index;
		regargs->distoparam.index = model->index;
		regargs->distoparam.filename = g_strdup(model->filename);
		regargs->ext_ref_undistorted = TRUE;
	}

	siril_log_message(_("LSF: registering %d frames against the Stack of %s\n"), seq->selnum, f->path);
	int ret = register_star_alignment(regargs);
	*layer = regargs->layer;
	free_regargs(regargs);
	if (ret || !seq->regparam[*layer]) {
		siril_log_error(_("LSF: registration against the Stack failed\n"));
		return 1;
	}
	return 0;
}

/* second pass: the footprint of the transforms of seq->regparam[layer]
 * decides the expansion of the grid (pad), and the registration is applied
 * onto the grid, expanded or not. Returns the registered sequence, written
 * with the lsf_r_ prefix. The transforms in seq stay relative to the grid
 * before the expansion. */
static sequence *apply_onto_grid(sequence *seq, const lsf_file *f, const struct lsf_cmd_args *a, int layer,
		const disto_params *model, int pad[4]) {
	const lsf_settings *s = &f->settings;
	if (!seq->regparam || !seq->regparam[layer]) {
		siril_log_error(_("LSF: the sequence has no registration on layer %d\n"), layer);
		return NULL;
	}
	/* Siril's apply undistorts the frames with the model recorded on the
	 * registration layer: the file's on an undistorted grid, none otherwise,
	 * whatever an earlier registration of the sequence left there (in memory
	 * only, the .seq file is not written) */
	if (!seq->distoparam && !(seq->distoparam = calloc(seq->nb_layers, sizeof(disto_params)))) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	seq->distoparam[layer].index = model ? model->index : DISTO_UNDEF;
	seq->distoparam[layer].filename = model ? g_strdup(model->filename) : NULL;

	/* footprint of the registered frames on the grid, without the edge
	 * samples the merge drops: the outermost one lies outside the frame's
	 * outer pixel centers already, each further one moves the edge in a pixel */
	double box[4];
	int first = -1;
	const double inset = MAX(lsf_edge_margin(s->interpolation) - 1, 0);
	lsf_footprint_init(box);
	for (int i = 0; i < seq->number; i++) {
		const Homography *H = &seq->regparam[layer][i].H;
		if (!seq->imgparam[i].incl || guess_transform_from_H(*H) == NULL_TRANSFORMATION)
			continue;
		double h[9] = { H->h00, H->h01, H->h02, H->h10, H->h11, H->h12, H->h20, H->h21, H->h22 };
		lsf_footprint_add(h, seq->is_variable ? seq->imgparam[i].rx : seq->rx,
				seq->is_variable ? seq->imgparam[i].ry : seq->ry, inset, box);
		if (first < 0)
			first = i;
	}
	if (first < 0) {
		siril_log_error(_("LSF: no frame could be registered against the Stack\n"));
		return NULL;
	}
	if (plan_expansion(f, a, box, pad))
		return NULL;

	/* Apply checks the reference frame of the sequence, which plays no part
	 * with an external reference: use a registered one. */
	struct registration_args *regargs = new_regargs(seq, s);
	if (!regargs)
		return NULL;
	regargs->func = register_apply_reg;
	regargs->layer = layer;
	regargs->no_output = FALSE;
	regargs->ext_ref_undistorted = model != NULL;
	regargs->keep_transforms = model && model->index == DISTO_FILES;	// -wcsreg, onto the Stack
	if (any_pad(pad))
		regargs->canvas =(framing_roi) { -pad[LSF_PAD_LEFT], -pad[LSF_PAD_LAST],
			f->rx + pad[LSF_PAD_LEFT] + pad[LSF_PAD_RIGHT], f->ry + pad[LSF_PAD_FIRST] + pad[LSF_PAD_LAST] };
	int seq_ref = seq->reference_image;
	seq->reference_image = first;
	regargs->reference_image = first;
	int ret = register_apply_reg(regargs);
	seq->reference_image = seq_ref;
	free_regargs(regargs);
	if (ret) {
		siril_log_error(_("LSF: applying the registration onto the grid failed\n"));
		return NULL;
	}
	gchar *rname = g_strdup_printf("%s%s", LSF_REG_PREFIX, seq->seqname);
	sequence *rseq = load_batch(rname);
	g_free(rname);
	return rseq;
}

/* both passes */
static sequence *register_against_stack(sequence *seq, const lsf_file *f, const struct lsf_cmd_args *a,
		const disto_params *model, int pad[4]) {
	int layer;
	if (compute_against_stack(seq, f, model, &layer))
		return NULL;
	return apply_onto_grid(seq, f, a, layer, model, pad);
}

/* first pass through the WCS (lsf_append -wcsreg, Section 10.4): each
 * frame's transform from its WCS and the Stack's. On an undistorted grid the
 * linear part of the frames' WCS is used, since the apply undistorts them. */
static int compute_against_stack_wcs(sequence *seq, const lsf_file *f, const struct lsf_cmd_args *a, int *layer) {
	const lsf_settings *s = &f->settings;
	if (!has_wcs((fits *) &f->primary)) {
		siril_log_error(_("LSF: the Stack of %s has no WCS, so frames cannot be registered through theirs. Plate-solve "
					"the Stack (lsf_export, then platesolve) and import its WCS with lsf_util -wcs=, or register "
					"with stars\n"), f->path);
		return 1;
	}
	*layer = seq->nb_layers == 3 ? 1 : 0;
	struct wcsprm *ref = wcs_copy_linear(f->primary.keywords.wcslib);
	if (!ref)
		return 1;
	if (f->primary.keywords.wcslib->lin.dispre && !s->undistorted) {
		/* the Stack keeps its own distortion: map through it */
		wcsfree(ref);
		free(ref);
		ref = wcs_deepcopy(f->primary.keywords.wcslib, NULL);
	}
	siril_log_message(_("LSF: registering %d frames against the Stack of %s through their WCS\n"), seq->selnum, f->path);
	/* the star measures are made on the grid, when the frames are merged
	 * (Section 10.3) */
	int ret = lsf_register_by_wcs(seq, *layer, ref, f->ry, s->transform, s->undistorted, a->max_threads);
	wcsfree(ref);
	free(ref);
	if (ret) {
		siril_log_error(_("LSF: no frame could be registered through its WCS: are the frames plate-solved?\n"));
		return ret;
	}
	/* the transforms map onto the Stack, as an external reference: Siril's
	 * apply takes them as they are, and gives the frames the Stack's WCS */
	seq->ext_ref = TRUE;
	g_free(seq->ext_ref_path);
	seq->ext_ref_path = g_strdup(f->path);
	seq->ext_ref_rx = f->rx;
	seq->ext_ref_ry = f->ry;
	return 0;
}

/***** registration made beforehand against an LSF file *****/

/* The binding: a sidecar file beside the .seq file, naming the LSF file a
 * sequence was registered against and its grid then. The .seq file holds
 * the transforms and the path of the external reference, but not what an
 * append needs to trust them. */
#define LSF_BINDING_GROUP "lsf_registration"

static gchar *binding_path(const sequence *seq) {
	return g_strdup_printf("%s.lsfreg", seq->seqname);
}

/* the transforms of seq on layer at the precision the .seq file stores them
 * (%g), so the fingerprint is the same before and after the file is written
 * and read: it changes when the sequence is registered again. Formatted
 * without the locale, so that a GUI and a CLI run agree. */
static gchar *transforms_fingerprint(const sequence *seq, int layer) {
	if (layer < 0 || layer >= seq->nb_layers || !seq->regparam || !seq->regparam[layer])
		return NULL;
	GChecksum *ck = g_checksum_new(G_CHECKSUM_SHA256);
	gchar buf[G_ASCII_DTOSTR_BUF_SIZE];
	for (int i = 0; i < seq->number; i++) {
		const Homography *H = &seq->regparam[layer][i].H;
		const double h[9] = { H->h00, H->h01, H->h02, H->h10, H->h11, H->h12, H->h20, H->h21, H->h22 };
		for (int k = 0; k < 9; k++) {
			g_ascii_formatd(buf, sizeof buf, "%g", h[k]);
			g_checksum_update(ck, (const guchar *) buf, strlen(buf));
			g_checksum_update(ck, (const guchar *) (k < 8 ? " " : "\n"), 1);
		}
	}
	gchar *hex = g_strdup(g_checksum_get_string(ck));
	g_checksum_free(ck);
	return hex;
}

/* what a registration against an LSF file keeps until its args are freed:
 * the hook_data of the registration */
struct lsf_reg_data {
	gchar *path, *uid;
	int gen, rx, ry, origin_x, origin_y;
	star_finder_params sf;	// the star detection stored in the file, regargs->sf_params
	gchar *disto;		// the distortion model the frames were undistorted with (JSON), or NULL
	gboolean registered;	// the registration succeeded
};

static void free_reg_data(gpointer p) {
	struct lsf_reg_data *d = p;
	if (!d)
		return;
	g_free(d->path);
	g_free(d->uid);
	g_free(d->disto);
	g_free(d);
}

/* the Living Stack the args registered against, if the registration
 * succeeded, or NULL; for the GUI, which then sets the stacking tab to
 * append to it. Valid until free_registration_hook_data(). */
const gchar *lsf_registration_result(const struct registration_args *regargs) {
	if (regargs->hook_free != free_reg_data || !regargs->hook_data)
		return NULL;
	const struct lsf_reg_data *d = regargs->hook_data;
	return d->registered ? d->path : NULL;
}

/* end_hook of the registration, in its thread after the .seq file is written */
static void lsf_registration_end(struct registration_args *regargs) {
	struct lsf_reg_data *d = regargs->hook_free == free_reg_data ? regargs->hook_data : NULL;
	if (!d)
		return;
	d->registered = !regargs->retval;
	gchar *fp = regargs->retval ? NULL : transforms_fingerprint(regargs->seq, regargs->layer);
	if (fp) {
		GKeyFile *kf = g_key_file_new();
		gchar *date = lsf_now_utc_iso();
		g_key_file_set_string(kf, LSF_BINDING_GROUP, "file", d->path);
		g_key_file_set_string(kf, LSF_BINDING_GROUP, "uid", d->uid);
		g_key_file_set_integer(kf, LSF_BINDING_GROUP, "generation", d->gen);
		g_key_file_set_integer(kf, LSF_BINDING_GROUP, "width", d->rx);
		g_key_file_set_integer(kf, LSF_BINDING_GROUP, "height", d->ry);
		g_key_file_set_integer(kf, LSF_BINDING_GROUP, "origin_x", d->origin_x);
		g_key_file_set_integer(kf, LSF_BINDING_GROUP, "origin_y", d->origin_y);
		g_key_file_set_integer(kf, LSF_BINDING_GROUP, "layer", regargs->layer);
		g_key_file_set_string(kf, LSF_BINDING_GROUP, "transforms", fp);
		if (d->disto)
			g_key_file_set_string(kf, LSF_BINDING_GROUP, "disto", d->disto);
		g_key_file_set_string(kf, LSF_BINDING_GROUP, "date", date);
		gchar *bpath = binding_path(regargs->seq);
		GError *err = NULL;
		if (g_key_file_save_to_file(kf, bpath, &err))
			siril_log_info(_("LSF: %s is registered against %s, generation %d. lsf_append uses this registration "
						"while the grid of the file does not change; do not apply it before appending\n"),
					regargs->seq->seqname, d->path, d->gen);
		else {
			siril_log_error(_("LSF: could not save %s: %s. lsf_append will register the frames again\n"),
					bpath, err ? err->message : "");
			g_clear_error(&err);
		}
		g_free(bpath);
		g_free(date);
		g_key_file_free(kf);
	}
	g_free(fp);
}

/* Prepares a registration (register command or Registration tab) whose
 * external reference is an LSF file: global star alignment, compute only,
 * no star lists, no distortion correction, with the transform, the minimum
 * pairs and the star detection stored in the file unless the caller set
 * them (keep_*). Its end writes the binding. Returns 0, or 1 with the
 * reason logged. */
int lsf_prepare_registration(struct registration_args *regargs, gboolean keep_type, gboolean keep_pairs,
		gboolean keep_findstar) {
	lsf_info info;
	const char *path = regargs->external_ref_path;
	if (lsf_read_info(path, &info)) {
		siril_log_error(_("LSF: %s cannot be registered against: %s\n"), path, info.problem);
		lsf_info_free(&info);
		return 1;
	}
	int ret = 1;
	const lsf_settings *s = &info.file.settings;
	if (info.no_append)
		siril_log_error(_("LSF: frames cannot be appended to %s with this version (%s), so registering against it "
					"is not useful\n"), path, info.no_append);
	else if (regargs->func != register_star_alignment)
		siril_log_error(_("LSF: only the global star alignment registers against a Living Stack, without the "
					"two-pass option\n"));
	else if (regargs->undistort != DISTO_UNDEF && !s->undistorted)
		siril_log_error(_("LSF: %s was created without distortion correction, and its frames cannot be undistorted: "
					"they would not line up with the Stack toward the edges. Register without -disto\n"), path);
	else if (s->undistorted && regargs->undistort == DISTO_IMAGE)
		siril_log_error(_("LSF: the distortion model of the current image cannot be used against a Living Stack, "
					"use -disto=file or -disto=master, or the model stored in the file\n"));
	else if (s->undistorted && regargs->undistort == DISTO_UNDEF && s->disto.index == DISTO_UNDEF)
		siril_log_error(_("LSF: %s is undistorted and stores no distortion model: give one with -disto=\n"), path);
	else if (s->undistorted && regargs->undistort == DISTO_UNDEF && s->disto.index == DISTO_FILES)
		siril_log_error(_("LSF: %s undistorts each frame with the SIP terms of its own WCS, which star alignment "
					"cannot use: plate-solve the frames and append them with lsf_append -wcsreg, or register "
					"with -disto=file or -disto=master\n"), path);
	else if (regargs->seq->nb_layers != info.file.nchan)
		siril_log_error(_("LSF: the Living Stack has %d channel(s) and the sequence %d\n"),
				info.file.nchan, regargs->seq->nb_layers);
	else ret = 0;
	if (ret) {
		lsf_info_free(&info);
		return ret;
	}
	struct lsf_reg_data *d = g_new0(struct lsf_reg_data, 1);
	d->path = g_canonicalize_filename(path, NULL);
	d->uid = g_strdup(info.file.uid);
	d->gen = info.file.gen;
	d->rx = info.file.rx;
	d->ry = info.file.ry;
	d->origin_x = s->origin_x;
	d->origin_y = s->origin_y;
	/* applying belongs to lsf_append, which grows the grid as needed */
	regargs->no_output = TRUE;
	regargs->no_starlist = TRUE;
	if (!keep_type)
		regargs->type = s->transform;
	if (!keep_pairs)
		regargs->min_pairs = s->min_pairs;
	if (!keep_findstar) {
		d->sf = com.pref.starfinder_conf;
		lsf_findstar_from_json(s->findstar_json, &d->sf);
		regargs->sf_params = &d->sf;
	}
	/* an undistorted grid: the frames are undistorted with the model given,
	 * or the file's, and the Stack, undistorted already, is not */
	if (s->undistorted) {
		if (regargs->undistort == DISTO_UNDEF) {
			regargs->undistort = s->disto.index;
			g_free(regargs->distoparam.filename);
			regargs->distoparam.index = s->disto.index;
			regargs->distoparam.filename = g_strdup(s->disto.filename);
		}
		if (regargs->distoparam.index == DISTO_MASTER && !regargs->distoparam.filename)
			regargs->distoparam.filename = g_strdup(com.pref.prepro.disto_lib);
		if (regargs->distoparam.index == DISTO_FILE && regargs->distoparam.filename) {
			gchar *abs = g_canonicalize_filename(regargs->distoparam.filename, NULL);	// as lsf_append stores it
			g_free(regargs->distoparam.filename);
			regargs->distoparam.filename = abs;
		}
		regargs->ext_ref_undistorted = TRUE;
		d->disto = lsf_disto_to_json(&regargs->distoparam);
		gchar *model = lsf_disto_describe(&regargs->distoparam);
		siril_log_message(_("LSF: the Living Stack is undistorted: the frames are undistorted with %s\n"), model);
		g_free(model);
	}
	regargs->hook_data = d;
	regargs->hook_free = free_reg_data;
	regargs->end_hook = lsf_registration_end;
	siril_log_message(_("LSF: registering against the Living Stack %s (generation %d, %dx%d grid): transforms "
				"only, applied by lsf_append\n"), path, info.file.gen, info.file.rx, info.file.ry);
	lsf_info_free(&info);
	return 0;
}

/* Whether the transforms of seq were computed against this file on its
 * current grid: the same LSFUID, grid size and origin, and
 * transforms unchanged since. A later generation with the same grid keeps
 * them valid; a grid grown since does not. Returns TRUE with the layer and,
 * if wanted, the date of the registration; otherwise FALSE with why, or no
 * why when the sequence was never registered against an LSF file. */
gboolean lsf_binding_check(const sequence *seq, const lsf_file *f, const disto_params *model, int *layer,
		gchar **date, gchar **why) {
	if (date)
		*date = NULL;
	if (why)
		*why = NULL;
	gchar *bpath = binding_path(seq);
	GKeyFile *kf = g_key_file_new();
	gboolean loaded = g_key_file_load_from_file(kf, bpath, G_KEY_FILE_NONE, NULL);
	g_free(bpath);
	if (!loaded) {
		g_key_file_free(kf);
		return FALSE;
	}
	gchar *uid = g_key_file_get_string(kf, LSF_BINDING_GROUP, "uid", NULL);
	gchar *fp = g_key_file_get_string(kf, LSF_BINDING_GROUP, "transforms", NULL);
	int rx = g_key_file_get_integer(kf, LSF_BINDING_GROUP, "width", NULL);
	int ry = g_key_file_get_integer(kf, LSF_BINDING_GROUP, "height", NULL);
	int ox = g_key_file_get_integer(kf, LSF_BINDING_GROUP, "origin_x", NULL);
	int oy = g_key_file_get_integer(kf, LSF_BINDING_GROUP, "origin_y", NULL);
	int l = g_key_file_get_integer(kf, LSF_BINDING_GROUP, "layer", NULL);
	gchar *bdisto = g_key_file_get_string(kf, LSF_BINDING_GROUP, "disto", NULL);
	gchar *mdisto = model && model->index != DISTO_UNDEF ? lsf_disto_to_json(model) : NULL;
	gchar *current = transforms_fingerprint(seq, l);
	const char *reason = NULL;
	if (!uid || strcmp(uid, f->uid))
		reason = _("the sequence was registered against another Living Stack");
	else if (rx != f->rx || ry != f->ry || ox != f->settings.origin_x || oy != f->settings.origin_y)
		reason = _("the grid of the Living Stack has grown since the sequence was registered against it");
	else if (!fp || !current || strcmp(fp, current))
		reason = _("the sequence was registered again since it was registered against the Living Stack");
	else if (g_strcmp0(bdisto, mdisto))
		reason = _("the sequence was registered with another distortion model than the append uses");
	gboolean ok = reason == NULL;
	if (ok) {
		*layer = l;
		if (date)
			*date = g_key_file_get_string(kf, LSF_BINDING_GROUP, "date", NULL);
	} else if (why)
		*why = g_strdup(reason);
	g_free(uid);
	g_free(fp);
	g_free(bdisto);
	g_free(mdisto);
	g_free(current);
	g_key_file_free(kf);
	return ok;
}

/* grows the grid of the file in memory before merging (Section 10.5): the
 * accumulator planes, grid.origin_*, the WCS of the Stack, and the shift
 * applied to the transforms of the ledger when the file is written */
static int expand_grid(lsf_file *f, const int pad[4], int gen) {
	int orx = f->rx, ory = f->ry;
	if (lsf_history_grow(&f->hist, pad))
		return 1;
	f->rx = f->hist.rx;
	f->ry = f->hist.ry;
	f->settings.origin_x += pad[LSF_PAD_LEFT];
	f->settings.origin_y += pad[LSF_PAD_FIRST];
	f->xfrm_shift[0] = pad[LSF_PAD_LEFT];
	f->xfrm_shift[1] = pad[LSF_PAD_LAST];
	/* the reference stars of the wFWHM stay on the same sky: in FITS pixels,
	 * as CRPIX */
	for (int i = 0; i < f->wref.n; i++) {
		f->wref.x[i] += pad[LSF_PAD_LEFT];
		f->wref.y[i] += pad[LSF_PAD_FIRST];
	}
	if (has_wcs(&f->primary)) {
		/* CRPIX counts pixels in storage order */
		Homography T;
		cvGetEye(&T);
		T.h02 = pad[LSF_PAD_LEFT];
		T.h12 = pad[LSF_PAD_FIRST];
		reframe_wcs(f->primary.keywords.wcslib, &T);
	}
	f->primary.history = g_slist_append(f->primary.history, g_strdup_printf(
				"LSF gen %d: grid %dx%d expanded to %dx%d, pad %d %d %d %d", gen, orx, ory,
				f->rx, f->ry, pad[LSF_PAD_LEFT], pad[LSF_PAD_RIGHT], pad[LSF_PAD_FIRST], pad[LSF_PAD_LAST]));
	return 0;
}

/* the registered frames carry the WCS that Siril's apply derived for the
 * expanded grid; it must agree with the one computed for the Stack */
static void check_expanded_wcs(sequence *rseq, fits *primary) {
	fits meta = { 0 };
	if (!has_wcs(primary) || seq_read_frame_metadata(rseq, 0, &meta))
		return;
	if (has_wcs(&meta)) {
		const double *a = primary->keywords.wcslib->crpix, *b = meta.keywords.wcslib->crpix;
		if (fabs(a[0] - b[0]) > 1e-6 || fabs(a[1] - b[1]) > 1e-6)
			siril_log_warning(_("LSF: the WCS of the expanded Stack (CRPIX %g %g) and of the registered "
						"frames (CRPIX %g %g) disagree\n"), a[0], a[1], b[0], b[1]);
		else siril_log_message(_("LSF: the WCS of the expanded Stack agrees with the registered frames\n"));
	}
	clearfits(&meta);
}

/* index in seq of each frame of the registered sequence rseq, or NULL */
static int *map_to_sources(sequence *seq, sequence *rseq) {
	int *src_index = malloc(rseq->number * sizeof(int));
	if (!src_index) {
		PRINT_ALLOC_ERR;
		return NULL;
	}
	for (int i = 0; i < rseq->number; i++) {
		src_index[i] = -1;
		for (int j = 0; j < seq->number; j++)
			if (seq->imgparam[j].filenum == rseq->imgparam[i].filenum) {
				src_index[i] = j;
				break;
			}
		if (src_index[i] < 0) {
			siril_log_error(_("LSF: cannot map registered frame %d back to its source\n"), rseq->imgparam[i].filenum);
			free(src_index);
			return NULL;
		}
	}
	return src_index;
}

/* lsf_init on a sequence whose registration is computed but not applied
 * (Section 10.1): applies it with Siril's apply step and the framing asked
 * for: onto the union of all the selected frames (max, the default; every
 * frame keeps the size of the grid instead of Siril's per-frame crop), the
 * reference frame (current), their intersection (min) or the reference
 * frame's size centred on their centre of gravity (cog). ref is the
 * reference frame, replaced by the first registered one when it is not
 * selected, origin the margins of the grid before it in storage order
 * (negative when the grid starts inside it) and base the transform from the
 * registration of the sequence onto the grid. Returns the registered
 * sequence, written with the lsf_r_ prefix. */
/* An astrometric registration (seqplatesolve): the transforms of the frames
 * from their WCS, by Siril's functions, for the framing asked. Max framing is
 * Siril's: a projection at the centre of the frames, at the angle that frames
 * them best, whichever way the reference frame is turned. The apply keeps
 * these transforms (keep_transforms) and gives the frames wcsref. Returns 0,
 * or 1 with the reason logged. */
static int astrometric_transforms(sequence *seq, const struct lsf_cmd_args *a, int layer, int ref,
		struct wcsprm **wcsref) {
	struct registration_args *regargs = new_regargs(seq, &a->settings);
	int *included = calloc(seq->number, sizeof(int));
	int ret = 1;
	if (!regargs || !included
			|| !(regargs->WCSDATA = calloc(seq->number, sizeof(struct wcsprm)))) {
		PRINT_ALLOC_ERR;
		goto end;
	}
	regargs->layer = layer;
	if (convert_parsed_filter_to_filter(&regargs->filters, seq, &regargs->filtering_criterion,
				&regargs->filtering_parameter) || collect_sequence_astrometry(regargs, included))
		goto end;
	if (!included[ref]) {
		siril_log_error(_("LSF: the reference frame of %s is not plate-solved\n"), seq->seqname);
		goto end;
	}
	int seq_ref = seq->reference_image;
	ret = compute_Hs_from_astrometry(seq, included, ref, regargs->WCSDATA, a->framing, layer, NULL, wcsref)
		|| !*wcsref;
	seq->reference_image = seq_ref;
end:
	if (ret)
		siril_log_error(_("LSF: the transforms of %s could not be computed from the WCS of its frames\n"), seq->seqname);
	if (regargs && regargs->WCSDATA && included) {
		for (int i = 0; i < seq->number; i++)
			if (included[i])
				wcsfree(regargs->WCSDATA + i);
		free(regargs->WCSDATA);
	}
	free(included);
	if (regargs)
		free_regargs(regargs);
	return ret;
}

static sequence *apply_init_registration(sequence *seq, const struct lsf_cmd_args *a, int *ref,
		int origin[2], Homography *base, gboolean *undistorted) {
	int layer = get_registration_layer(seq);
	int pad[4] = { 0 };
	if (seq->ext_ref) {
		siril_log_error(_("LSF: the registration of %s was computed against an external reference image. "
					"lsf_init can only apply a registration onto one of the frames: register the sequence "
					"again, or apply its registration first\n"), seq->seqname);
		return NULL;
	}
	gboolean ok = seq->imgparam[*ref].incl && guess_transform_from_H(seq->regparam[layer][*ref].H) != NULL_TRANSFORMATION;
	for (int i = 0; !ok && i < seq->number; i++) {
		if (seq->imgparam[i].incl && guess_transform_from_H(seq->regparam[layer][i].H) != NULL_TRANSFORMATION) {
			siril_log_warning(_("LSF: the reference image is not in the selected set of images, "
						"using image %d as temporary reference image\n"), i + 1);
			*ref = i;
			ok = TRUE;
		}
	}
	if (!ok) {
		siril_log_error(_("LSF: no selected frame of %s has registration data\n"), seq->seqname);
		return NULL;
	}
	gboolean astrometric = layer_has_distortion(seq, layer) && seq->distoparam[layer].index == DISTO_FILES;
	struct wcsprm *wcsref = NULL;
	if (astrometric && astrometric_transforms(seq, a, layer, *ref, &wcsref))
		return NULL;
	Homography Href = seq->regparam[layer][*ref].H;
	/* what the apply maps the frames by: Href^-1 * H, and H as it is for an
	 * astrometric registration, whose transforms are onto the grid's WCS */
	Homography Hframe = Href;
	if (astrometric)
		cvGetEye(&Hframe);
	int rx = seq->is_variable ? seq->imgparam[*ref].rx : seq->rx;
	int ry = seq->is_variable ? seq->imgparam[*ref].ry : seq->ry;

	if (a->framing == FRAMING_MAX) {
		/* footprint of the frames on the reference frame, as a canvas */
		double box[4];
		lsf_footprint_init(box);
		for (int i = 0; i < seq->number; i++) {
			Homography Hi = seq->regparam[layer][i].H, Hr;
			if (!seq->imgparam[i].incl || guess_transform_from_H(Hi) == NULL_TRANSFORMATION)
				continue;
			cvTransfH(&Hi, &Hframe, &Hr);	// Hframe^-1 * Hi
			double h[9] = { Hr.h00, Hr.h01, Hr.h02, Hr.h10, Hr.h11, Hr.h12, Hr.h20, Hr.h21, Hr.h22 };
			lsf_footprint_add(h, seq->is_variable ? seq->imgparam[i].rx : seq->rx,
					seq->is_variable ? seq->imgparam[i].ry : seq->ry, 0.0, box);
		}
		lsf_grid_overhang(box, rx, ry, pad);
	}
	struct registration_args *regargs = new_regargs(seq, &a->settings);
	if (!regargs)
		return NULL;
	regargs->func = register_apply_reg;
	regargs->layer = layer;
	regargs->no_output = FALSE;
	regargs->framing = a->framing == FRAMING_MAX ? FRAMING_CURRENT : a->framing;
	regargs->keep_transforms = astrometric;
	regargs->wcsref = wcsref;
	if (any_pad(pad))
		regargs->canvas = (framing_roi) { -pad[LSF_PAD_LEFT], -pad[LSF_PAD_LAST],
			rx + pad[LSF_PAD_LEFT] + pad[LSF_PAD_RIGHT], ry + pad[LSF_PAD_FIRST] + pad[LSF_PAD_LAST] };
	int seq_ref = seq->reference_image;
	seq->reference_image = *ref;
	regargs->reference_image = *ref;
	int ret = register_apply_reg(regargs);
	seq->reference_image = seq_ref;
	if (ret) {
		free_regargs(regargs);
		siril_log_error(_("LSF: applying the registration of %s failed\n"), seq->seqname);
		return NULL;
	}
	/* Siril's apply maps each frame onto the grid by Htransf^-1 * H, for
	 * every framing; the reference frame's first stored pixel gives origin */
	Homography R, Ht = regargs->framingd.Htransf;
	int nrx = regargs->framingd.roi_out.w, nry = regargs->framingd.roi_out.h;
	*undistorted = regargs->undistort != DISTO_UNDEF;
	free_regargs(regargs);
	*base = Ht;
	cvInvertH(base);
	cvMultH(*base, Href, &R);
	const double hr[9] = { R.h00, R.h01, R.h02, R.h10, R.h11, R.h12, R.h20, R.h21, R.h22 };
	lsf_frame_origin(hr, rx, ry, nry, origin);
	siril_log_message(_("LSF: registration of %s applied with %s framing: grid %dx%d, origin %d %d "
				"from the reference frame\n"), seq->seqname, lsf_framing_name(a->framing), nrx, nry,
			origin[0], origin[1]);
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
	/* -out writes the appended stack to a new file and leaves the source as
	 * it was, as lsf_util -out does */
	const char *target = a->out ? a->out : a->path;
	if (a->out && g_file_test(a->out, G_FILE_TEST_EXISTS) && !a->overwrite) {
		siril_log_error(_("LSF: %s already exists, use -overwrite to replace it\n"), a->out);
		return 1;
	}
	GPtrArray *rows = g_ptr_array_new_with_free_func(lsf_ledger_row_free);

	if (lsf_open_and_validate(a->path, &f, FALSE)) {
		lsf_file_free(&f);
		g_ptr_array_free(rows, TRUE);
		return 1;
	}
	if (a->out)
		f.primary.history = g_slist_append(f.primary.history, g_strdup_printf(
					"LSF generation %d: written from %s", f.gen + 1, f.path));
	/* this writer merges onto a grid at the frames' own scale only: frames
	 * merged onto a drizzled grid (Section 6) would land on the wrong pixels */
	if (f.settings.oversample != 1) {
		siril_log_error(_("LSF: %s has a grid drizzled at scale %d, which this version cannot append to. "
					"Only grid.oversample 1 is supported\n"), a->path, f.settings.oversample);
		goto end;
	}
	/* overrides of the defaults stored in the file apply from this append on
	 * and are saved for the following ones (Section 6). A HISTORY card
	 * records each change with the generation that made it. */
	double hist_before[3] = { f.settings.hist_low, f.settings.hist_high, f.settings.hist_margin };
	apply_hist_overrides(a, &f.settings);
	if (hist_before[0] != f.settings.hist_low || hist_before[1] != f.settings.hist_high
			|| hist_before[2] != f.settings.hist_margin) {
		siril_log_message(_("LSF: history test changed from low %g high %g margin %g to low %g high %g margin %g. "
					"It applies from this append on and is saved in the file for the following ones; "
					"the frames already in the stack are not affected\n"),
				hist_before[0], hist_before[1], hist_before[2],
				f.settings.hist_low, f.settings.hist_high, f.settings.hist_margin);
		f.primary.history = g_slist_append(f.primary.history, g_strdup_printf(
					"LSF gen %d: history test low %g high %g margin %g, was %g %g %g", f.gen + 1,
					f.settings.hist_low, f.settings.hist_high, f.settings.hist_margin,
					hist_before[0], hist_before[1], hist_before[2]));
	}
	/* without rejection the two values mean nothing */
	if (a->rej_set && (a->settings.rej != f.settings.rej || (a->settings.rej != NO_REJEC
					&& (a->settings.sig[0] != f.settings.sig[0] || a->settings.sig[1] != f.settings.sig[1])))) {
		gchar *was = rejection_text(f.settings.rej, f.settings.sig), *now = rejection_text(a->settings.rej, a->settings.sig);
		siril_log_message(_("LSF: rejection changed from %s to %s. It applies from this append on and is saved in "
					"the file for the following ones; the frames already in the stack are not affected\n"), was, now);
		f.primary.history = g_slist_append(f.primary.history, g_strdup_printf(
					"LSF gen %d: rejection %s, was %s", f.gen + 1, now, was));
		g_free(was);
		g_free(now);
		f.settings.rej = a->settings.rej;
		f.settings.sig[0] = a->settings.sig[0];
		f.settings.sig[1] = a->settings.sig[1];
	}
	if (a->rejmaps && f.settings.rejmaps)
		siril_log_message(_("LSF: %s already has rejection maps, counted since generation %d\n"),
				a->path, f.hist.rej_gen);
	f.settings.rejmaps |= a->rejmaps;	// lsf_process_batch() creates them
	/* distortion correction is a property of the grid (Section 10.4): a grid
	 * created with it undistorts every batch, with the model stored in the
	 * file or the one given, which becomes the new default; a grid created
	 * without it never does */
	const disto_params *model = NULL;
	if (a->disto.index != DISTO_UNDEF && !f.settings.undistorted) {
		siril_log_error(_("LSF: %s was created without distortion correction and cannot switch: its Stack keeps the "
					"distortion of its first batch, and undistorted frames would not line up with it toward the "
					"edges. Append without -disto\n"), a->path);
		goto end;
	}
	if (f.settings.undistorted && !a->noreg) {
		if (a->disto.index != DISTO_UNDEF && !lsf_disto_same(&a->disto, &f.settings.disto)) {
			gchar *was = lsf_disto_describe(&f.settings.disto), *now = lsf_disto_describe(&a->disto);
			siril_log_message(_("LSF: distortion model changed from %s to %s. It applies from this append on and is "
						"saved in the file for the following ones\n"), was, now);
			f.primary.history = g_slist_append(f.primary.history, g_strdup_printf(
						"LSF gen %d: distortion model %s, was %s", f.gen + 1, now, was));
			g_free(was);
			g_free(now);
			lsf_disto_copy(&f.settings.disto, &a->disto);
		}
		if (f.settings.disto.index == DISTO_UNDEF) {
			siril_log_error(_("LSF: %s is undistorted and stores no distortion model: give one with -disto=\n"), a->path);
			goto end;
		}
		/* the frames' own SIP terms: Siril's star alignment undistorts the
		 * stars with one model, so only the WCS registers such frames */
		if (f.settings.disto.index == DISTO_FILES && !a->wcsreg) {
			siril_log_error(_("LSF: %s undistorts each frame with the SIP terms of its own WCS (it was created from an "
						"astrometric registration): plate-solve the frames and append with -wcsreg, or give one "
						"model for all of them with -disto=file or -disto=master\n"), a->path);
			goto end;
		}
		model = &f.settings.disto;
	} else if (f.settings.undistorted)
		siril_log_warning(_("LSF: %s is undistorted: -noreg merges the frames as they are, which must be undistorted "
					"already, or they will not line up with the Stack toward the edges\n"), a->path);
	/* the weight model switches from range to inverse from this append on
	 * (Section 10.3): the merged frames keep their weights */
	if (a->weight_inverse) {
		lsf_settings *s = &f.settings;
		if (s->weighting != WFWHM_WEIGHT && s->weighting != NBSTARS_WEIGHT) {
			siril_log_error(_("LSF: -weight_model applies to wfwhm and nbstars weighting, which freeze a range; %s is "
						"weighted by %s\n"), a->path, lsf_weight_name(s->weighting));
			goto end;
		}
		if (s->weight_model == LSF_WMODEL_INVERSE)
			siril_log_message(_("LSF: %s already uses the inverse weight model\n"), a->path);
		else {
			float *hist_wfwhm = NULL;
			int *hist_nstars = NULL;
			double *hist_weight = NULL;
			long n = lsf_read_ledger_weights(&f, &hist_wfwhm, &hist_nstars, &hist_weight);
			int failed = n <= 0 || lsf_switch_to_inverse(s, (int) n, hist_wfwhm, hist_nstars, hist_weight);
			free(hist_wfwhm);
			free(hist_nstars);
			free(hist_weight);
			if (failed) {
				siril_log_error(_("LSF: the ledger of %s does not give the star measures of its frames, the weight "
							"model cannot be switched\n"), a->path);
				goto end;
			}
			siril_log_message(_("LSF: weight model switched from range to inverse, from this append on: frames beyond "
						"the range of the first batch are no longer clamped. The %ld frames already merged keep "
						"their weights, and the new ones are on the same scale on average\n"), n);
			f.primary.history = g_slist_append(f.primary.history, g_strdup_printf(
						"LSF gen %d: weight model inverse, was range", f.gen + 1));
		}
	}
	seq = load_batch(a->seqname);
	if (!seq)
		goto end;
	if (seq->nb_layers != f.nchan) {
		siril_log_error(_("LSF: the file has %d channel(s), the batch has %d\n"), f.nchan, seq->nb_layers);
		goto end;
	}
	gboolean listed = apply_frame_list(seq, a);

	sequence *batch = seq;
	int gen = f.gen + 1, batchnum = f.nb_batches + 1;
	int pad[4] = { 0 };
	if (!a->noreg) {
		int layer = 0;
		if (a->wcsreg) {
			/* through the WCS: a registration made beforehand is not used */
			if (!compute_against_stack_wcs(seq, &f, a, &layer))
				rseq = apply_onto_grid(seq, &f, a, layer, model, pad);
		} else {
			/* a registration made beforehand against this file
			 * replaces the first pass, unless -reg asks for a new one */
			gchar *date = NULL, *why = NULL;
			gboolean reuse = !a->force_reg && lsf_binding_check(seq, &f, model, &layer, &date, &why);
			if (reuse) {
				siril_log_message(_("LSF: using the registration of %s against this file, made %s\n"),
						a->seqname, date ? date : "");
				rseq = apply_onto_grid(seq, &f, a, layer, model, pad);
			} else {
				if (why && !a->force_reg)
					siril_log_message(_("LSF: %s: registering the frames again\n"), why);
				rseq = register_against_stack(seq, &f, a, model, pad);
			}
			g_free(date);
			g_free(why);
		}
		if (!rseq)
			goto end;
		if (any_pad(pad)) {
			if (expand_grid(&f, pad, gen))
				goto end;
			check_expanded_wcs(rseq, &f.primary);
		}
		if (!(src_index = map_to_sources(seq, rseq)))
			goto end;
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

	struct lsf_batch_args bargs = { 0 };
	bargs.seq = batch;
	bargs.src_seq = a->noreg ? NULL : seq;
	bargs.src_index = src_index;
	bargs.framing = a->noreg ? NULL : lsf_framing_name(a->framing);
	bargs.reg_method = a->wcsreg ? "wcs" : "stars";
	bargs.disto = model;
	bargs.settings = &f.settings;
	bargs.hist = &f.hist;
	bargs.init = FALSE;
	bargs.batch = batchnum;
	bargs.gen = gen;
	bargs.prev_frames = f.ledger_rows;
	bargs.prev_exptime = f.ledger_exptime;
	bargs.seq_name = seq->seqname;
	bargs.rows = rows;
	bargs.primary = &f.primary;
	bargs.wref = &f.wref;
	bargs.keep_seq_file = listed && batch == seq;
	bargs.max_threads = a->max_threads;
	bargs.overlap_norm = a->overlap_norm;
	bargs.feather = a->feather;
	if (lsf_process_batch(&bargs))
		goto end;
	/* the transforms of the new rows were computed onto the grid before its expansion */
	for (guint i = 0; i < rows->len && (f.xfrm_shift[0] != 0.0 || f.xfrm_shift[1] != 0.0); i++) {
		lsf_ledger_row *r = g_ptr_array_index(rows, i);
		lsf_xfrm_shift(r->xfrm, f.xfrm_shift[0], f.xfrm_shift[1]);
	}
	if (lsf_build_primary(&f.primary, &f.settings, &f.hist))
		goto end;
	/* the file keeps its compression unless the command changes it */
	lsf_compression comp = a->compression_set ? a->compression : f.compression;
	if (comp != f.compression)
		siril_log_message(_("LSF: compression of the extensions changed from %s to %s\n"),
				lsf_compression_name(f.compression), lsf_compression_name(comp));
	gchar *replace_msg = NULL;
	int written = lsf_write_file(target, &f, &f.settings, &f.hist, &f.wref, &f.primary, f.uid, gen, rows, batchnum,
			a->backup && !a->out, comp, &replace_msg);
	if (note_write(a, written, replace_msg))
		goto end;
	if (a->out)
		siril_log_info(_("LSF: %u frames appended to %s and written to %s, %ld frames in total, generation %d. "
					"%s is unchanged: append to one of the two files only\n"), rows->len, a->path, a->out,
				f.ledger_rows + rows->len, gen, a->path);
	else siril_log_info(_("LSF: %u frames appended to %s, %ld frames in total, generation %d\n"),
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

static void print_summary(const lsf_file *f) {
	const lsf_settings *s = &f->settings;
	const fkeywords *k = &f->primary.keywords;
	siril_log_message(_("LSF file %s: version %s, id %s, generation %d\n"), f->path, LSF_VERSION, f->uid, f->gen);
	siril_log_message(_("  %ld frames in %d batch(es), STACKCNT %u, live time %.1f s\n"), f->ledger_rows,
			f->nb_batches, k->stackcnt, k->livetime > 0.0 ? k->livetime : 0.0);
	if (s->is_cfa)
		siril_log_message(_("  grid %dx%d, %d channel(s), CFA sensor %s (debayered)\n"), f->rx, f->ry, f->nchan, s->bayer);
	else siril_log_message(_("  grid %dx%d, %d channel(s), mono sensor\n"), f->rx, f->ry, f->nchan);
	if (has_wcs((fits *) &f->primary)) {
		gchar *desc = lsf_wcs_describe(f->primary.keywords.wcslib, f->rx, f->ry);
		siril_log_message(_("  WCS: %s\n"), desc);
		g_free(desc);
	} else siril_log_message(_("  no WCS\n"));
	if (s->undistorted) {
		gchar *model = lsf_disto_describe(&s->disto);
		siril_log_message(_("  undistorted grid, distortion model %s\n"), model);
		g_free(model);
	}
	gchar *rej = rejection_text(s->rej, s->sig);
	siril_log_message(_("  rejection %s, history test low %g high %g margin %g\n"),
			rej, s->hist_low, s->hist_high, s->hist_margin);
	g_free(rej);
	siril_log_message(_("  normalization %s%s%s, output normalization %s\n"), lsf_norm_name(s->norm),
			s->lite_norm ? _(" (fast)") : "", s->rgb_equal ? _(" (RGB equalized)") : "", s->outnorm ? _("on") : _("off"));
	switch (s->weighting) {
	case NOISE_WEIGHT:
		siril_log_message(_("  weighting noise, frozen norm %g %g %g\n"), s->noise_norm[0],
				f->nchan == 3 ? s->noise_norm[1] : 0.0, f->nchan == 3 ? s->noise_norm[2] : 0.0);
		break;
	case WFWHM_WEIGHT:
		siril_log_message(_("  weighting wfwhm (%s), frozen min %g max %g norm %g, floor %g, range x%g\n"),
				lsf_weight_model_name(s->weight_model), s->fwhm_min, s->fwhm_max, s->fwhm_norm, s->weight_floor,
				s->weight_range);
		break;
	case NBSTARS_WEIGHT:
		siril_log_message(_("  weighting nbstars (%s), frozen min %g max %g norm %g, floor %g, range x%g\n"),
				lsf_weight_model_name(s->weight_model), s->stars_min, s->stars_max, s->stars_norm, s->weight_floor,
				s->weight_range);
		break;
	default:
		siril_log_message(_("  weighting none\n"));
	}
	double raw = (double) f->rx * f->ry * f->nchan;
	siril_log_message(_("  LSF_ACC %s, %.0f MB stored (%.0f%% of raw), LSF_COUNT %s, %.0f MB stored (%.0f%% of raw)\n"),
			f->acc_codec, (double) f->acc_bytes / BYTES_IN_A_MB,
			100.0 * f->acc_bytes / (raw * (LSF_NB_PL + (f->has_exp ? 1 : 0)) * sizeof(float)),
			f->count_codec, (double) f->count_bytes / BYTES_IN_A_MB, 100.0 * f->count_bytes / (raw * sizeof(WORD)));
	if (f->has_exp && f->hist.exp[0]) {
		float emax = 0.f;
		for (int c = 0; c < f->nchan; c++)
			for (size_t i = 0; i < (size_t) f->rx * f->ry; i++)
				emax = MAX(emax, f->hist.exp[c][i]);
		siril_log_message(_("  exposure map (EXP planes of LSF_ACC), up to %.1f s per pixel\n"), emax);
	} else if (f->has_exp)
		siril_log_message(_("  exposure map (EXP planes of LSF_ACC)\n"));
	else if (isfinite(f->ledger_exptime))
		siril_log_message(_("  no exposure map: every frame has %g s, so COUNT x %g s is the exposure per pixel\n"),
				f->ledger_exptime, f->ledger_exptime);
	else siril_log_message(_("  no exposure map, and the exposures of the frames are not all the same or not all known: "
				"the exposure per pixel is not known\n"));
	if (s->rejmaps)	// LSF_REJ is present, check 7
		siril_log_message(_("  LSF_REJ %s, %.0f MB stored (%.0f%% of raw), rejections counted since generation %d%s\n"),
				f->rej_codec, (double) f->rej_bytes / BYTES_IN_A_MB, 100.0 * f->rej_bytes / (raw * 2 * sizeof(WORD)),
				f->hist.rej_gen, f->hist.rej_saturated ? ", some counts saturated" : "");
	else siril_log_message(_("  no rejection maps (LSF_REJ)\n"));
}

/* lsf_util: Section 8 validation and a summary (-check), Repair or Adopt
 * when only check 8 failed, and a change of the compression of the
 * extensions (Section 5.3). The actions combine into one rewrite, which
 * increments the generation. */
/* lsf_util -info: the summary from the headers alone (lsf_read_info) */
static int lsf_info_run(struct lsf_cmd_args *a) {
	lsf_info info;
	if (lsf_read_info(a->path, &info)) {
		siril_log_error(_("LSF: %s cannot be used, %s\n"), a->path, info.problem);
		lsf_info_free(&info);
		return 1;
	}
	print_summary(&info.file);
	siril_log_message(_("  %.1f s of exposure, frames from %s to %s, last commit %s, file %.0f MB\n"), info.exposure,
			info.first_date ? info.first_date : "?", info.last_date ? info.last_date : "?",
			info.last_added ? info.last_added : "?", (double) info.file_bytes / BYTES_IN_A_MB);
	if (info.no_append)
		siril_log_warning(_("LSF: this version cannot append to %s: %s\n"), a->path, info.no_append);
	siril_log_info(_("LSF: %s passed checks 1 to 4, 6 and 7 of Section 8; -check runs them all\n"), a->path);
	lsf_info_free(&info);
	return 0;
}

/* The checks lsf_append makes of a batch, for the GUI before it runs:
 * a channel count that differs is refused, and frames whose Section 7
 * keywords differ from the file are skipped, which the reference frame
 * shows. Returns a message for the user or NULL, *fatal when refused. */
gchar *lsf_check_batch(const lsf_info *info, sequence *seq, gboolean *fatal) {
	*fatal = FALSE;
	if (seq->nb_layers != info->file.nchan) {
		*fatal = TRUE;
		return g_strdup_printf(_("The Living Stack has %d channel(s) and the sequence %d."),
				info->file.nchan, seq->nb_layers);
	}
	fits meta = { 0 };
	gchar *msg = NULL;
	if (!seq_read_frame_metadata(seq, sequence_find_refimage(seq), &meta)) {
		const char *bad = lsf_frame_mismatch(&info->file.primary, &meta);
		if (bad)
			msg = g_strdup_printf(_("The %s of the reference frame differs from the Living Stack. Frames "
						"that differ are not merged."), bad);
	}
	clearfits(&meta);
	return msg;
}

/* free space of the file system holding the directory path, as Siril
 * measures it (find_space), or -1. fs_id identifies that file system, so that
 * two directories on the same one are not counted twice. */
static gint64 free_space_of(const char *path, gchar **fs_id) {
	if (fs_id) {
		GFile *file = g_file_new_for_path(path);
		GFileInfo *fi = g_file_query_info(file, G_FILE_ATTRIBUTE_ID_FILESYSTEM, G_FILE_QUERY_INFO_NONE, NULL, NULL);
		if (fi) {
			*fs_id = g_strdup(g_file_info_get_attribute_string(fi, G_FILE_ATTRIBUTE_ID_FILESYSTEM));
			g_object_unref(fi);
		}
		g_object_unref(file);
	}
	return find_space(path);
}

/* Before an append, from lsf_read_info: memory for the accumulator, which
 * stays in memory during the merge, and disk for the registered frames
 * (written in the working directory, unless noreg) and for the new file,
 * written beside the old one before it replaces it (Section 10.2), or as
 * out when not NULL. The grid may still grow, which is known only after the
 * registration. Returns a message or NULL, *fatal when the append cannot
 * run. */
gchar *lsf_append_preflight(const lsf_info *info, const char *wd, const char *out, int nb_frames, gboolean noreg,
		gboolean add_rejmaps, gboolean *fatal) {
	const lsf_file *f = &info->file;
	*fatal = FALSE;
	gboolean rej = f->settings.rejmaps || add_rejmaps;
	guint64 npix = (guint64) f->rx * f->ry * f->nchan;
	/* an exposure map the append creates is not counted: the batch is not read yet */
	guint64 hist = lsf_history_bytes(npix, rej, f->has_exp);
	guint64 max_mem = (guint64) get_max_memory_in_MB() * BYTES_IN_A_MB;
	if (hist >= max_mem) {
		*fatal = TRUE;
		return g_strdup_printf(_("The accumulator of this Living Stack needs %.0f MB of memory, and Siril may use "
					"%.0f MB. Raise the memory limit in the preferences."),
				(double) hist / BYTES_IN_A_MB, (double) max_mem / BYTES_IN_A_MB);
	}
	GString *msg = g_string_new(NULL);
	if (hist > max_mem / 4 * 3)
		g_string_append_printf(msg, _("The accumulator needs %.0f MB of the %.0f MB Siril may use, which leaves "
					"little memory for the frames: the append may be slow or fail. "),
				(double) hist / BYTES_IN_A_MB, (double) max_mem / BYTES_IN_A_MB);

	/* the new file is about the size of the old one, plus the growth */
	gint64 need_file = info->file_bytes;
	gint64 need_frames = noreg ? 0 : (gint64) nb_frames * (gint64) (npix * sizeof(float));
	gchar *dir = g_path_get_dirname(out ? out : f->path), *fs_file = NULL, *fs_wd = NULL;
	gint64 free_file = free_space_of(dir, &fs_file), free_wd = free_space_of(wd, &fs_wd);
	g_free(dir);
	gboolean same_fs = fs_file && fs_wd && !strcmp(fs_file, fs_wd);
	gboolean short_file = free_file >= 0 && free_file < need_file + (same_fs ? need_frames : 0);
	gboolean short_wd = !same_fs && free_wd >= 0 && free_wd < need_frames;
	if (short_file || short_wd) {
		*fatal = TRUE;
		gchar *nf = g_format_size(need_file), *nr = g_format_size(need_frames);
		if (out)
			g_string_printf(msg, _("Not enough free disk space: the append writes the registered frames (%s, in the "
						"working directory) and the new Living Stack file (%s)."), nr, nf);
		else g_string_printf(msg, _("Not enough free disk space: the append writes the registered frames (%s, in the "
					"working directory) and a new copy of the Living Stack (%s, beside it) before replacing it."), nr, nf);
		g_free(nf);
		g_free(nr);
	}
	g_free(fs_file);
	g_free(fs_wd);
	return msg->len ? g_string_free(msg, FALSE) : (g_string_free(msg, TRUE), NULL);
}

/* Before lsf_init, from the sequence alone: memory for the accumulator of
 * the new grid, and disk for the registered frames, written in the working
 * directory when the registration is not applied yet, and for the file.
 * The grid is known only once the registration is applied: the size of the
 * frames stands for it, and max framing makes it larger. Returns a message
 * or NULL, *fatal when the creation cannot run. */
gchar *lsf_init_preflight(sequence *seq, const char *wd, const char *path, int nb_frames, gboolean rejmaps,
		gboolean *fatal) {
	*fatal = FALSE;
	int nchan = seq->nb_layers;
	guint64 npix = (guint64) seq->rx * seq->ry * nchan;
	guint64 hist = lsf_history_bytes(npix, rejmaps, FALSE);
	guint64 max_mem = (guint64) get_max_memory_in_MB() * BYTES_IN_A_MB;
	if (hist >= max_mem) {
		*fatal = TRUE;
		return g_strdup_printf(_("The accumulator of this Living Stack needs at least %.0f MB of memory, and Siril may "
					"use %.0f MB. Raise the memory limit in the preferences."),
				(double) hist / BYTES_IN_A_MB, (double) max_mem / BYTES_IN_A_MB);
	}
	GString *msg = g_string_new(NULL);
	if (hist > max_mem / 4 * 3)
		g_string_append_printf(msg, _("The accumulator needs at least %.0f MB of the %.0f MB Siril may use, which "
					"leaves little memory for the frames: the creation may be slow or fail. "),
				(double) hist / BYTES_IN_A_MB, (double) max_mem / BYTES_IN_A_MB);

	/* the file: the Stack (float32) and the accumulator, uncompressed */
	gint64 need_file = (gint64) (npix * sizeof(float) + hist);
	gboolean applied = !layer_has_usable_registration(seq, get_registration_layer(seq));
	gint64 need_frames = applied ? 0 : (gint64) nb_frames * (gint64) (npix * sizeof(float));
	gchar *dir = g_path_get_dirname(path), *fs_file = NULL, *fs_wd = NULL;
	gint64 free_file = free_space_of(dir, &fs_file), free_wd = free_space_of(wd, &fs_wd);
	g_free(dir);
	gboolean same_fs = fs_file && fs_wd && !strcmp(fs_file, fs_wd);
	gboolean short_file = free_file >= 0 && free_file < need_file + (same_fs ? need_frames : 0);
	gboolean short_wd = !same_fs && free_wd >= 0 && free_wd < need_frames;
	if (short_file || short_wd) {
		*fatal = TRUE;
		gchar *nf = g_format_size(need_file), *nr = g_format_size(need_frames);
		g_string_printf(msg, _("Not enough free disk space: the creation writes the registered frames (%s, in the "
					"working directory) and the Living Stack (about %s) at least."), nr, nf);
		g_free(nf);
		g_free(nr);
	}
	g_free(fs_file);
	g_free(fs_wd);
	return msg->len ? g_string_free(msg, FALSE) : (g_string_free(msg, TRUE), NULL);
}

static int lsf_util_run(struct lsf_cmd_args *a) {
	lsf_file f;
	int retval = 1;
	GPtrArray *none = NULL;
	fits imported = { 0 };
	struct wcsprm *new_wcs = NULL;
	if (a->info)
		return lsf_info_run(a);
	gchar *problem = NULL;
	if (lsf_open_and_validate_msg(a->path, &f, TRUE, &problem)) {
		add_message(a, _("The file failed the checks: %s. It was not modified."), problem ? problem : _("see the log"));
		g_free(problem);
		lsf_file_free(&f);
		return 1;
	}
	print_summary(&f);
	a->stack_modified = f.stack_modified;
	if (!f.stack_modified) {
		siril_log_info(_("LSF: %s passed all the checks of Section 8\n"), a->path);
		if (a->check)
			add_message(a, _("The file passed all the checks."));
	}
	gboolean fix_stack = (a->repair || a->adopt) && f.stack_modified;
	if ((a->repair || a->adopt) && !f.stack_modified) {
		siril_log_message(_("LSF: the Stack matches the accumulator, there is nothing to repair or adopt\n"));
		add_message(a, _("The Stack matches the accumulator, there is nothing to repair or adopt."));
	}
	gboolean recompress = a->compression_set && a->compression != f.compression;
	if (a->compression_set && !recompress) {
		if (f.compression == LSF_COMP_NONE)
			siril_log_message(_("LSF: the extensions of %s are already uncompressed\n"), a->path);
		else siril_log_message(_("LSF: the extensions of %s are already compressed with %s\n"), a->path,
				lsf_compression_name(f.compression));
		add_message(a, _("The compression is already %s."), lsf_compression_name(f.compression));
	}
	/* any rewrite rebuilds the Stack from MEAN, which would silently repair it */
	if (f.stack_modified && !fix_stack) {
		siril_log_error(_("LSF: %s fails check 8, its Stack was modified outside Siril. lsf_append will "
					"refuse it, and lsf_util will not rewrite it, until lsf_util is run with -repair (rebuild "
					"the Stack from the accumulator) or -adopt (keep the edit)\n"), a->path);
		add_message(a, _("The Stack (the image of the file) was modified outside Siril; the accumulator is intact. "
					"Appending is refused until the Stack is repaired, which undoes the edit, or adopted, "
					"which keeps it."));
		retval = a->check && !a->repair && !a->adopt && !a->compression_set && !a->wcs_file && !a->out ? 0 : 1;
		goto end;
	}
	/* WCS import (Section 8): the WCS of HDU 0, which the writer copies to
	 * every extension, is replaced; the pixels and LSFHASH do not change */
	if (a->wcs_file) {
		gchar *why = NULL;
		if (lsf_wcs_read_import(a->wcs_file, &f, &imported, &why)) {
			siril_log_error(_("LSF: the WCS cannot be imported: %s\n"), why);
			add_message(a, _("The WCS cannot be imported: %s."), why);
			g_free(why);
			goto end;
		}
		new_wcs = wcs_deepcopy(imported.keywords.wcslib, NULL);
		/* an undistorted grid keeps a linear WCS: SIP terms solved on it only
		 * fit residuals, and would make Siril distort the Stack's stars */
		if (new_wcs && f.settings.undistorted && new_wcs->lin.dispre) {
			double shift = lsf_wcs_sip_shift(new_wcs, f.rx, f.ry);
			struct wcsprm *lin = wcs_copy_linear(new_wcs);
			wcsfree(new_wcs);
			free(new_wcs);
			new_wcs = lin;
			siril_log_message(_("LSF: the grid is undistorted, so the SIP terms of the imported WCS are dropped. They "
						"corrected at most %.2f px at the edges of the grid\n"), shift);
			if (shift > 0.5)
				siril_log_warning(_("LSF: a correction of %.2f px on an undistorted grid suggests that the distortion "
							"model does not fit these frames well\n"), shift);
		}
		if (!new_wcs) {
			siril_log_error(_("LSF: the WCS of %s cannot be copied\n"), a->wcs_file);
			goto end;
		}
		if (has_wcs(&f.primary) && lsf_wcs_same(new_wcs, f.primary.keywords.wcslib, f.rx, f.ry)) {
			siril_log_message(_("LSF: %s already has this WCS\n"), a->path);
			add_message(a, _("The file already has this WCS."));
			wcsfree(new_wcs);
			free(new_wcs);
			new_wcs = NULL;
		}
	}
	/* -out writes the result to a new file and leaves the source untouched,
	 * even when no action changes anything */
	const char *target = a->out ? a->out : a->path;
	if (!fix_stack && !recompress && !new_wcs && !a->out) {
		retval = 0;
		goto end;
	}
	if (a->out && g_file_test(a->out, G_FILE_TEST_EXISTS) && !a->overwrite) {
		siril_log_error(_("LSF: %s already exists, use -overwrite to replace it\n"), a->out);
		goto end;
	}
	gchar *card;
	int gen = f.gen + 1;
	size_t adopted = 0, adopted_nonfinite = 0;
	if (a->out)
		f.primary.history = g_slist_append(f.primary.history, g_strdup_printf(
					"LSF generation %d: written from %s", gen, f.path));
	if (fix_stack) {
		if (a->adopt) {
			size_t n = (size_t) f.rx * f.ry, changed = 0, nonfinite = 0;
			gboolean norm = f.settings.outnorm;
			if (norm && !(f.settings.outnorm_max > f.settings.outnorm_min)) {
				siril_log_error(_("LSF: the output normalization range of %s is not valid, the edited Stack cannot be adopted\n"), a->path);
				goto end;
			}
			for (int c = 0; c < f.nchan; c++)
				changed += lsf_adopt_plane(f.hist.plane[LSF_PL_MEAN][c], f.hist.plane[LSF_PL_W][c],
						f.stack_data + c * n, n, norm, f.settings.outnorm_min, f.settings.outnorm_max, &nonfinite);
			siril_log_message(_("LSF: %zu edited pixel(s) adopted into MEAN, M2 is unchanged\n"), changed);
			if (nonfinite)
				siril_log_warning(_("LSF: %zu edited pixel(s) are NaN or infinite and were not adopted, "
							"they keep the value of the accumulator\n"), nonfinite);
			card = g_strdup_printf("LSF generation %d: adopted the Stack edited outside Siril, %zu pixels", gen, changed);
			adopted = changed;
			adopted_nonfinite = nonfinite;
		} else {
			card = g_strdup_printf("LSF generation %d: Stack edited outside Siril, repaired from MEAN", gen);
		}
		f.primary.history = g_slist_append(f.primary.history, card);
	}
	lsf_compression comp = a->compression_set ? a->compression : f.compression;
	if (recompress) {
		if (comp == LSF_COMP_NONE)
			card = g_strdup_printf("LSF generation %d: extensions decompressed", gen);
		else card = g_strdup_printf("LSF generation %d: extensions compressed with %s", gen, lsf_compression_name(comp));
		f.primary.history = g_slist_append(f.primary.history, card);
	}
	gboolean wcs_imported = new_wcs != NULL;
	if (new_wcs) {
		gchar *desc = lsf_wcs_describe(new_wcs, f.rx, f.ry), *base = g_path_get_basename(a->wcs_file);
		f.primary.history = g_slist_append(f.primary.history, g_strdup_printf(
					"LSF generation %d: WCS imported from %s", gen, base));
		siril_log_message(_("LSF: WCS imported from %s: %s\n"), a->wcs_file, desc);
		g_free(desc);
		g_free(base);
		free_wcs(&f.primary);
		f.primary.keywords.wcslib = new_wcs;
		new_wcs = NULL;
		/* the centre keywords Siril writes with a WCS */
		f.primary.keywords.wcsdata = imported.keywords.wcsdata;
		center2wcs2(f.primary.keywords.wcslib, f.rx, f.ry, &f.primary.keywords.wcsdata.ra,
				&f.primary.keywords.wcsdata.dec);
	}
	if (lsf_build_primary(&f.primary, &f.settings, &f.hist))
		goto end;
	none = g_ptr_array_new();
	gchar *replace_msg = NULL;
	int written = lsf_write_file(target, &f, &f.settings, &f.hist, &f.wref, &f.primary, f.uid, gen, none, f.nb_batches,
			a->backup, comp, &replace_msg);
	if (note_write(a, written, replace_msg))
		goto end;
	if (fix_stack && !a->out)
		a->stack_modified = FALSE;	// the file now passes check 8
	if (fix_stack && a->adopt) {
		siril_log_info(_("LSF: edited Stack of %s adopted, generation %d\n"), target, gen);
		add_message(a, _("The edited Stack was adopted: %zu pixels changed the mean, the variance is unchanged."),
				adopted);
		if (adopted_nonfinite)
			add_message(a, _("%zu edited pixels are NaN or infinite and kept the value of the accumulator."),
					adopted_nonfinite);
	} else if (fix_stack) {
		siril_log_info(_("LSF: Stack of %s repaired from the accumulator, generation %d\n"), target, gen);
		add_message(a, _("The Stack was rebuilt from the accumulator."));
	}
	if (recompress && comp == LSF_COMP_NONE) {
		siril_log_info(_("LSF: extensions of %s decompressed, generation %d\n"), target, gen);
		add_message(a, _("The extensions are no longer compressed."));
	} else if (recompress) {
		siril_log_info(_("LSF: extensions of %s compressed with %s, generation %d\n"), target,
				lsf_compression_name(comp), gen);
		add_message(a, _("The extensions are compressed with %s."), lsf_compression_name(comp));
	}
	if (wcs_imported) {
		siril_log_info(_("LSF: WCS of %s imported from %s, generation %d\n"), target, a->wcs_file, gen);
		add_message(a, _("The WCS was written into the Living Stack."));
	}
	if (a->out) {
		siril_log_info(_("LSF: %s written from %s, generation %d, %s left unchanged\n"), a->out, a->path,
				gen, a->path);
		add_message(a, _("Written to %s; the original file is unchanged."), a->out);
	}
	retval = 0;
end:
	if (none)
		g_ptr_array_free(none, TRUE);
	if (new_wcs) {
		wcsfree(new_wcs);
		free(new_wcs);
	}
	clearfits(&imported);
	lsf_file_free(&f);
	return retval;
}

/* lsf_maps: images derived from the accumulator (Sections 10.7, 10.8). The
 * noise and sigma maps are in the units of the Stack, the exposure map in
 * seconds. */
static int lsf_maps_run(struct lsf_cmd_args *a) {
	static const struct { guint flag; int rej; const char *suffix, *what; } kinds[] = {
		{ LSF_MAP_NEFF, -1, "neff", N_("effective number of frames") },
		{ LSF_MAP_NOISE, -1, "noise", N_("standard deviation of the mean") },
		{ LSF_MAP_SIGMA, -1, "sigma", N_("standard deviation of one frame") },
		{ LSF_MAP_SNR, -1, "snr", N_("mean over its standard deviation") },
		{ LSF_MAP_REJLO, 0, "rejlo", N_("samples rejected low") },
		{ LSF_MAP_REJHI, 1, "rejhi", N_("samples rejected high") },
		{ LSF_MAP_EXP, -1, "exp", N_("exposure of the merged samples, seconds") },
	};
	lsf_file f;
	int retval = 0;
	if (lsf_open_and_validate(a->path, &f, TRUE)) {
		lsf_file_free(&f);
		return 1;
	}
	if (f.stack_modified)
		siril_log_warning(_("LSF: the maps come from the accumulator, which is intact\n"));
	if ((a->maps & (LSF_MAP_REJLO | LSF_MAP_REJHI)) && !f.hist.rej[0][0]) {
		siril_log_error(_("LSF: %s has no rejection maps. They are created by lsf_init or lsf_append "
					"with -rejmaps, and count from then on\n"), a->path);
		lsf_file_free(&f);
		return 1;
	}
	/* without the exposure map, all the frames had one exposure, or the
	 * exposure per pixel is not known (Section 10.8) */
	gboolean exp_derived = !f.hist.exp[0];
	if ((a->maps & LSF_MAP_EXP) && exp_derived && !isfinite(f.ledger_exptime)) {
		siril_log_error(_("LSF: %s has no exposure map, and the exposures of its frames are not all the same or not all "
					"known: the exposure per pixel cannot be derived\n"), a->path);
		lsf_file_free(&f);
		return 1;
	}
	const lsf_settings *s = &f.settings;
	double scale = s->outnorm && s->outnorm_max > s->outnorm_min ? 1.0 / (s->outnorm_max - s->outnorm_min) : 1.0;
	size_t n = (size_t) f.rx * f.ry;
	for (int m = 0; m < (int) G_N_ELEMENTS(kinds) && !retval; m++) {
		if (!(a->maps & kinds[m].flag))
			continue;
		fits *fit = NULL;
		int rej = kinds[m].rej;
		if (new_fit_image(&fit, f.rx, f.ry, f.nchan, rej >= 0 ? DATA_USHORT : DATA_FLOAT)) {
			retval = 1;
			break;
		}
		copy_fits_metadata(&f.primary, fit);
		lsf_strip_lsf_keys(fit);
		fit->keywords.lo = fit->keywords.hi = 0;
		fit->keywords.flo = fit->keywords.fhi = 0.f;
		/* rejection maps: the counts themselves, as LSF_REJ stores them */
		for (int c = 0; rej >= 0 && c < f.nchan; c++)
			memcpy(fit->pdata[c], f.hist.rej[rej][c], n * sizeof(WORD));
		gboolean exp = kinds[m].flag == LSF_MAP_EXP;
		for (int c = 0; exp && c < f.nchan; c++) {
			if (!exp_derived)
				memcpy(fit->fpdata[c], f.hist.exp[c], n * sizeof(float));
			else lsf_exp_backfill(fit->fpdata[c], f.hist.count[c], n, f.ledger_exptime);
		}
		for (int c = 0; rej < 0 && !exp && c < f.nchan; c++) {
			const float *mean = f.hist.plane[LSF_PL_MEAN][c], *W = f.hist.plane[LSF_PL_W][c];
			const float *M2 = f.hist.plane[LSF_PL_M2][c], *W2 = f.hist.plane[LSF_PL_W2][c];
			const WORD *N = f.hist.count[c];
			float *out = fit->fpdata[c];
			for (size_t i = 0; i < n; i++) {
				double neff, s2, vm, v;
				lsf_pixel_stats(s->weighting, W[i], M2[i], W2[i], N[i], &neff, &s2, &vm);
				switch (kinds[m].flag) {
				case LSF_MAP_NEFF: v = neff > 0.0 ? neff : 0.0; break;
				case LSF_MAP_NOISE: v = vm >= 0.0 ? sqrt(vm) * scale : 0.0; break;
				case LSF_MAP_SIGMA: v = s2 >= 0.0 ? sqrt(s2) * scale : 0.0; break;
				default: v = vm > 0.0 ? mean[i] / sqrt(vm) : 0.0; break;
				}
				out[i] = (float) v;
			}
		}
		if (rej >= 0)
			fit->history = g_slist_append(fit->history, g_strdup_printf("LSF %s map (%s since generation %d)",
						kinds[m].suffix, kinds[m].what, f.hist.rej_gen));
		else if (exp && exp_derived)
			fit->history = g_slist_append(fit->history, g_strdup_printf("LSF %s map (%s, COUNT x %g s)",
						kinds[m].suffix, kinds[m].what, f.ledger_exptime));
		else fit->history = g_slist_append(fit->history, g_strdup_printf("LSF %s map (%s)",
					kinds[m].suffix, kinds[m].what));
		fit->history = g_slist_append(fit->history, g_strdup_printf("from LSF file %s, generation %d",
					f.uid, f.gen));
		gchar *name = g_strdup_printf("%s_%s", a->out, kinds[m].suffix);
		if (savefits(name, fit)) {
			siril_log_error(_("LSF: could not save the %s map\n"), kinds[m].suffix);
			retval = 1;
		} else siril_log_message(_("LSF: %s map (%s) saved as %s\n"), kinds[m].suffix, _(kinds[m].what), name);
		g_free(name);
		clearfits(fit);
		free(fit);
	}
	lsf_file_free(&f);
	return retval;
}

static int lsf_export_run(const char *path, const char *out) {
	int retval = lsf_export_master(path, out);
	if (!retval) {
		/* the name savefits gave the file */
		gchar *name = set_right_extension(out);
		siril_log_message(_("LSF: master of %s exported to %s\n"), path, name ? name : out);
		g_free(name);
	}
	return retval;
}

/* the new Stack becomes the result image, as stacking leaves its result in
 * gfit with the sequence still loaded. It is a view of the LSF file, which
 * is never saved over (Section 13), or of the plain image exported with it. */
static void lsf_load_result(const char *path, const char *export_name) {
	fits fit = { 0 };
	if (readfits(export_name ? export_name : path, &fit, NULL, TRUE)) {
		clearfits(&fit);
		return;
	}
	lsf_strip_lsf_keys(&fit);
	single *uniq = calloc(1, sizeof(single));
	if (!uniq) {
		PRINT_ALLOC_ERR;
		clearfits(&fit);
		return;
	}
	g_rw_lock_writer_lock(&gfit->rwlock);
	clearfits(gfit);
	memcpy(gfit, &fit, offsetof(fits, rwlock));
	com.seq.current = RESULT_IMAGE;
	/* as stacking does, the previous com.uniq is not freed: the main thread
	 * may still be reading it */
	uniq->filename = strdup(export_name ? export_name : path);
	uniq->fileexist = TRUE;
	uniq->is_lsf = !export_name;
	uniq->comment = strdup(_("Living Stack result image"));
	uniq->nb_layers = gfit->naxes[2];
	uniq->fit = gfit;
	com.uniq = uniq;
	g_rw_lock_writer_unlock(&gfit->rwlock);
	/* with the writer lock released, as after stacking */
	notify_gfit_data_modified();
}

/* on the main thread when a job with an end function is over */
static gboolean lsf_end_job(gpointer p) {
	lsf_job_result *r = p;
	stop_processing_thread();
	r->end_func(r);
	g_free(r->path);
	g_free(r->out);
	g_free(r->export_out);
	g_free(r->message);
	free(r);
	return FALSE;
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
	case LSF_CMD_UTIL:
		retval = lsf_util_run(a);
		break;
	case LSF_CMD_MAPS:
		retval = lsf_maps_run(a);
		break;
	default:
		retval = lsf_export_run(a->path, a->out);
		break;
	}
	/* the plain image of the new Stack, as lsf_export makes it */
	gboolean commit = a->mode == LSF_CMD_INIT || a->mode == LSF_CMD_APPEND;
	const char *written = a->mode == LSF_CMD_APPEND && a->out ? a->out : a->path;
	gchar *export_name = commit && a->export_out ? set_right_extension(a->export_out) : NULL;
	if (!retval && export_name)
		retval = lsf_export_run(written, a->export_out);
	if (!retval && commit && a->load_result)
		lsf_load_result(written, export_name);
	lsf_job_result *r = a->end_func ? calloc(1, sizeof(lsf_job_result)) : NULL;
	if (r) {
		r->mode = a->mode;
		r->retval = retval;
		r->path = g_strdup(a->path);
		if (a->mode == LSF_CMD_APPEND || a->mode == LSF_CMD_UTIL)
			r->out = g_strdup(a->out);
		r->export_out = export_name;
		r->message = a->message;
		a->message = NULL;
		r->stack_modified = a->stack_modified;
		r->write_status = a->write_status;
		r->end_func = a->end_func;
		r->end_data = a->end_data;
	} else g_free(export_name);
	lsf_cmd_args_free(a);
	if (r)
		siril_add_idle(lsf_end_job, r);
	else siril_add_idle(end_generic, NULL);
	return GINT_TO_POINTER(retval);
}

/* starts the job of a in the processing thread, which takes it. Returns
 * non-zero, with a freed, when another job is running. */
int lsf_start(struct lsf_cmd_args *a) {
	if (!start_in_new_thread(lsf_worker, a)) {
		lsf_cmd_args_free(a);
		return 1;
	}
	return 0;
}

/* arg as parse_line() reads it back: quoted when it has blanks or starts with
 * a quote, with the quote it does not contain, as there is no escape */
static void add_arg(GString *cmd, const char *arg) {
	if (!strpbrk(arg, " \t") && arg[0] != '"' && arg[0] != '\'')
		g_string_append_printf(cmd, " %s", arg);
	else if (!strchr(arg, '"'))
		g_string_append_printf(cmd, " \"%s\"", arg);
	else g_string_append_printf(cmd, " '%s'", arg);
}

static void add_hist_and_common_args(GString *cmd, const struct lsf_cmd_args *a) {
	if (a->hist_low_set)
		g_string_append_printf(cmd, " -hist_low=%g", a->hist_low);
	if (a->hist_high_set)
		g_string_append_printf(cmd, " -hist_high=%g", a->hist_high);
	if (a->hist_margin_set)
		g_string_append_printf(cmd, " -hist_margin=%g", a->hist_margin);
	if (a->framing_set || a->framing != FRAMING_MAX)
		g_string_append_printf(cmd, " -framing=%s", lsf_framing_name(a->framing));
	if (a->rejmaps)
		g_string_append(cmd, " -rejmaps");
	if (a->feather > 0)
		g_string_append_printf(cmd, " -feather=%d", a->feather);
	if (a->compression_set)
		g_string_append(cmd, a->compression == LSF_COMP_NONE ? " -nocompress" :
				a->compression == LSF_COMP_GZIP1 ? " -compress=gzip1" : " -compress");
}

/* the command line of a, for the log of a GUI run, so it can be repeated
 * in a script. Returns NULL for the modes that have no GUI yet. */
gchar *lsf_cmd_args_to_command(const struct lsf_cmd_args *a) {
	GString *cmd = g_string_new(NULL);
	const lsf_settings *s = &a->settings;
	switch (a->mode) {
	case LSF_CMD_INIT:
		g_string_append(cmd, "lsf_init");
		add_arg(cmd, a->seqname);
		add_arg(cmd, a->path);
		if (s->rej != NO_REJEC)
			g_string_append_printf(cmd, " rej %s %g %g", rejection_type_letter(s->rej), s->sig[0], s->sig[1]);
		if (s->norm != NO_NORM)
			g_string_append_printf(cmd, " -norm=%s", lsf_norm_name(s->norm));
		if (s->lite_norm)
			g_string_append(cmd, " -fastnorm");
		if (s->rgb_equal)
			g_string_append(cmd, " -rgb_equal");
		if (s->weighting != NO_WEIGHT)
			g_string_append_printf(cmd, " -weight=%s", lsf_weight_name(s->weighting));
		if ((s->weighting == WFWHM_WEIGHT || s->weighting == NBSTARS_WEIGHT) && s->weight_floor != LSF_WEIGHT_FLOOR)
			g_string_append_printf(cmd, " -weight_floor=%g", s->weight_floor);
		if ((s->weighting == WFWHM_WEIGHT || s->weighting == NBSTARS_WEIGHT) && s->weight_range != 1.0)
			g_string_append_printf(cmd, " -weight_range=%g", s->weight_range);
		if (s->outnorm)
			g_string_append(cmd, " -output_norm");
		if (a->overlap_norm)
			g_string_append(cmd, " -overlap_norm");
		add_hist_and_common_args(cmd, a);
		if (a->overwrite)
			g_string_append(cmd, " -overwrite");
		break;
	case LSF_CMD_APPEND:
		g_string_append(cmd, "lsf_append");
		add_arg(cmd, a->path);
		add_arg(cmd, a->seqname);
		if (a->rej_set) {
			if (s->rej == NO_REJEC)
				g_string_append(cmd, " rej n");
			else g_string_append_printf(cmd, " rej %s %g %g", rejection_type_letter(s->rej), s->sig[0], s->sig[1]);
		}
		if (a->noreg)
			g_string_append(cmd, " -noreg");
		else if (a->wcsreg)
			g_string_append(cmd, " -wcsreg");
		if (a->weight_inverse)
			g_string_append(cmd, " -weight_model=inverse");
		else if (a->force_reg)
			g_string_append(cmd, " -reg");
		if (a->overlap_norm)
			g_string_append(cmd, " -overlap_norm");
		if (a->disto.index == DISTO_MASTER)
			g_string_append(cmd, " -disto=master");
		else if (a->disto.index == DISTO_FILES)
			g_string_append(cmd, " -disto=frames");
		else if (a->disto.index == DISTO_FILE) {
			g_string_append(cmd, " -disto=file");
			add_arg(cmd, a->disto.filename);
		}
		add_hist_and_common_args(cmd, a);
		if (a->expand_step != LSF_EXPAND_STEP)
			g_string_append_printf(cmd, " -expand_step=%d", a->expand_step);
		if (a->expand_force)
			g_string_append(cmd, " -expand_force");
		if (a->out) {
			gchar *opt = g_strdup_printf("-out=%s", a->out);
			add_arg(cmd, opt);
			g_free(opt);
			if (a->overwrite)
				g_string_append(cmd, " -overwrite");
		} else if (a->backup)
			g_string_append(cmd, " -backup");
		break;
	case LSF_CMD_EXPORT:
		g_string_append(cmd, "lsf_export");
		add_arg(cmd, a->path);
		add_arg(cmd, a->out);
		break;
	case LSF_CMD_UTIL:
		g_string_append(cmd, "lsf_util");
		add_arg(cmd, a->path);
		if (a->info)
			g_string_append(cmd, " -info");
		if (a->check)
			g_string_append(cmd, " -check");
		if (a->repair)
			g_string_append(cmd, " -repair");
		else if (a->adopt)
			g_string_append(cmd, " -adopt");
		if (a->compression_set)
			g_string_append(cmd, a->compression == LSF_COMP_NONE ? " -decompress" :
					a->compression == LSF_COMP_GZIP1 ? " -compress=gzip1" : " -compress");
		if (a->wcs_file) {
			gchar *opt = g_strdup_printf("-wcs=%s", a->wcs_file);
			add_arg(cmd, opt);
			g_free(opt);
		}
		if (a->out) {
			gchar *opt = g_strdup_printf("-out=%s", a->out);
			add_arg(cmd, opt);
			g_free(opt);
			if (a->overwrite)
				g_string_append(cmd, " -overwrite");
		} else if (a->backup)
			g_string_append(cmd, " -backup");
		break;
	default:
		g_string_free(cmd, TRUE);
		return NULL;
	}
	if (a->export_out && (a->mode == LSF_CMD_INIT || a->mode == LSF_CMD_APPEND)) {
		g_string_append(cmd, "\nlsf_export");
		add_arg(cmd, a->path);
		add_arg(cmd, a->export_out);
	}
	return g_string_free(cmd, FALSE);
}

/* -name=value with a value above 0, into out and set: 1 when consumed, 0 when
 * word_i is another option, -1 on a wrong value */
static int parse_double_opt(const char *word_i, const char *prefix, double *out, gboolean *set) {
	if (!g_str_has_prefix(word_i, prefix))
		return 0;
	const char *value = word_i + strlen(prefix);
	gchar *end;
	double v = g_ascii_strtod(value, &end);
	/* thresholds in sigma units: 0 would reject every sample on that side of the mean */
	if (end == value || *end || !(v > 0.0) || !isfinite(v)) {
		siril_log_error(_("Wrong parameter values. %s is not valid, aborting.\n"), word_i);
		return -1;
	}
	*out = v;
	*set = TRUE;
	return 1;
}

static int set_compression_opt(struct lsf_cmd_args *a, lsf_compression comp) {
	if (a->compression_set) {
		siril_log_error(_("LSF: give only one compression option\n"));
		return -1;
	}
	a->compression = comp;
	a->compression_set = TRUE;
	return 1;
}

/* -compress[=gzip2|gzip1] (Section 5.3): the accumulator is float data, which
 * only gzip compresses losslessly; returns 1 if consumed, -1 on error */
static int parse_compress_option(struct lsf_cmd_args *a, const char *w) {
	if (strcmp(w, "-compress") && !g_str_has_prefix(w, "-compress="))
		return 0;
	const char *v = w[9] == '=' ? w + 10 : "gzip2";
	if (!g_ascii_strcasecmp(v, "gzip2"))
		return set_compression_opt(a, LSF_COMP_GZIP2);
	if (!g_ascii_strcasecmp(v, "gzip1"))
		return set_compression_opt(a, LSF_COMP_GZIP1);
	siril_log_error(_("LSF: compression %s is not supported. The accumulator must be compressed losslessly, "
				"which only gzip1 and gzip2 do\n"), v);
	return -1;
}

/* options common to lsf_init and lsf_append; returns 1 if consumed, -1 on error */
static int parse_common_option(struct lsf_cmd_args *a, const char *w) {
	int c = parse_double_opt(w, "-hist_low=", &a->hist_low, &a->hist_low_set);
	if (!c)
		c = parse_double_opt(w, "-hist_high=", &a->hist_high, &a->hist_high_set);
	if (!c && (c = parse_double_opt(w, "-hist_margin=", &a->hist_margin, &a->hist_margin_set)) > 0
			&& a->hist_margin < 1.0) {
		siril_log_error(_("-hist_margin must be at least 1\n"));
		return -1;
	}
	if (c)
		return c;
	if (g_str_has_prefix(w, "-framing=")) {
		const char *v = w + 9;
		if (!strcmp(v, "max")) a->framing = FRAMING_MAX;
		else if (!strcmp(v, "current")) a->framing = FRAMING_CURRENT;
		else if (!strcmp(v, "min")) a->framing = FRAMING_MIN;
		else if (!strcmp(v, "cog")) a->framing = FRAMING_COG;
		else {
			siril_log_error(_("LSF: framing %s is not supported (max, current, min, cog)\n"), v);
			return -1;
		}
		a->framing_set = TRUE;
		return 1;
	}
	if (!strcmp(w, "-rejmaps")) {
		a->rejmaps = TRUE;
		return 1;
	}
	if (g_str_has_prefix(w, "-feather=")) {
		gchar *end;
		gint64 v = g_ascii_strtoll(w + 9, &end, 10);
		if (end == w + 9 || *end || v < 0 || v > LSF_FEATHER_MAX) {
			siril_log_error(_("Wrong parameter values. %s is not valid (0 to %d), aborting.\n"), w, LSF_FEATHER_MAX);
			return -1;
		}
		a->feather = (int) v;
		return 1;
	}
	if (!strcmp(w, "-nocompress"))
		return set_compression_opt(a, LSF_COMP_NONE);
	return parse_compress_option(a, w);
}

/* [rej [type] sigma_low sigma_high] at word[*pi], as the stack command takes
 * it, into the rejection of a->settings, advancing *pi. Returns 1 when
 * given, 0 when not, -1 on error. */
static int parse_rejection(struct lsf_cmd_args *a, int *pi) {
	lsf_settings *s = &a->settings;
	int i = *pi;
	if (!word[i] || (strcmp(word[i], "rej") && strcmp(word[i], "mean")))
		return 0;
	i++;
	int shift = 1;
	if (!parse_rejection_type(word[i], TRUE, &s->rej)) {
		s->rej = WINSORIZED;	// as the stack command, without a type
		shift = 0;
	}
	i += shift;
	if (s->rej != NO_REJEC || (word[i] && string_is_a_number(word[i]))) {
		gchar *end1 = NULL, *end2 = NULL;
		if (!word[i] || !word[i + 1] || !string_is_a_number(word[i]) || !string_is_a_number(word[i + 1])) {
			siril_log_error(_("The average stacking with rejection requires two extra arguments: sigma low and high.\n"));
			return -1;
		}
		s->sig[0] = g_ascii_strtod(word[i], &end1);
		s->sig[1] = g_ascii_strtod(word[i + 1], &end2);
		if (s->sig[0] < 0.0 || s->sig[1] < 0.0) {
			siril_log_error(_("The average stacking with rejection requires two extra arguments: sigma low and high.\n"));
			return -1;
		}
		if (!check_rejection_values(s->rej, s->sig[0], s->sig[1]))
			return -1;
		i += 2;
	}
	*pi = i;
	return 1;
}

/* lsf_init sequencename file.lsf.fits [rej [type] sigma_low sigma_high] [-nonorm, -norm=] [-fastnorm]
 *   [-rgb_equal] [-weight=] [-output_norm] [-framing=max|current|min|cog] [-hist_low=] [-hist_high=] [-hist_margin=]
 *   [-rejmaps] [-overlap_norm] [-feather=] [-compress[=gzip2|gzip1]] [-overwrite] */
int process_lsf_init(int nb) {
	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_INIT);
	if (!a)
		return CMD_ALLOC_ERROR;
	a->seqname = g_strdup(word[1]);
	a->path = g_strdup(word[2]);
	lsf_settings *s = &a->settings;

	int i = 3;
	if (parse_rejection(a, &i) < 0)
		goto failure;
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
		else if (!strcmp(w, "-overlap_norm")) a->overlap_norm = TRUE;
		else if (!strcmp(w, "-overwrite")) a->overwrite = TRUE;
		else if (g_str_has_prefix(w, "-weight_floor=")) {
			gchar *end;
			double v = g_ascii_strtod(w + 14, &end);
			if (end == w + 14 || *end || !(v > 0.0 && v < 1.0)) {
				siril_log_error(_("Wrong parameter values. %s is not valid (above 0 and below 1), aborting.\n"), w);
				goto failure;
			}
			s->weight_floor = v;
		}
		else if (g_str_has_prefix(w, "-weight_range=")) {
			gchar *end;
			double v = g_ascii_strtod(w + 14, &end);
			if (end == w + 14 || *end || !(v >= 1.0) || !isfinite(v)) {
				siril_log_error(_("Wrong parameter values. %s is not valid (at least 1), aborting.\n"), w);
				goto failure;
			}
			s->weight_range = v;
		}
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
	if (lsf_start(a))
		return CMD_GENERIC_ERROR;
	return CMD_OK;
failure:
	lsf_cmd_args_free(a);
	return CMD_ARG_ERROR;
}

/* the same file under two spellings, so that -out= naming the source is an
 * update in place */
gboolean lsf_same_file(const char *p1, const char *p2) {
	gchar *c1 = g_canonicalize_filename(p1, NULL), *c2 = g_canonicalize_filename(p2, NULL);
#ifdef _WIN32
	gchar *f1 = g_utf8_casefold(c1, -1), *f2 = g_utf8_casefold(c2, -1);
	gboolean same = !strcmp(f1, f2);
	g_free(f1);
	g_free(f2);
#else
	gboolean same = !strcmp(c1, c2);
#endif
	g_free(c1);
	g_free(c2);
	return same;
}

/* lsf_append file.lsf.fits sequencename [rej [type] sigma_low sigma_high] [-noreg] [-framing=max|current]
 *   [-noexpand] [-expand_step=] [-expand_force] [-out=file.lsf.fits [-overwrite] | -backup] [-hist_low=]
 *   [-hist_high=] [-hist_margin=] [-rejmaps] [-overlap_norm] [-feather=] [-compress[=gzip2|gzip1] | -nocompress] */
int process_lsf_append(int nb) {
	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_APPEND);
	if (!a)
		return CMD_ALLOC_ERROR;
	a->path = g_strdup(word[1]);
	a->seqname = g_strdup(word[2]);
	int i = 3;
	int rej = parse_rejection(a, &i);	// overrides the rejection stored in the file
	if (rej < 0)
		goto failure;
	a->rej_set = rej > 0;
	for (; i < nb && word[i]; i++) {
		const char *w = word[i];
		int c = parse_common_option(a, w);
		if (c < 0) goto failure;
		if (c > 0) continue;
		if (!strcmp(w, "-noreg")) a->noreg = TRUE;
		else if (!strcmp(w, "-reg")) a->force_reg = TRUE;
		else if (!strcmp(w, "-wcsreg")) a->wcsreg = TRUE;
		else if (!strcmp(w, "-weight_model=inverse")) a->weight_inverse = TRUE;
		else if (!strcmp(w, "-overlap_norm")) a->overlap_norm = TRUE;
		else if (!strcmp(w, "-disto=master")) {
			/* as register -disto=master: the pattern of the preferences */
			if (!com.pref.prepro.disto_lib || !com.pref.prepro.disto_lib[0]) {
				siril_log_error(_("-disto=master needs a distortion master pattern in the preferences\n"));
				goto failure;
			}
			g_free(a->disto.filename);
			a->disto.index = DISTO_MASTER;
			a->disto.filename = g_strdup(com.pref.prepro.disto_lib);
		}
		else if (!strcmp(w, "-disto=frames")) {
			/* the SIP terms of each frame's own WCS, as an astrometric
			 * registration undistorts them (-wcsreg) */
			g_free(a->disto.filename);
			a->disto = (disto_params) { .index = DISTO_FILES };
		}
		else if (!strcmp(w, "-disto=file")) {
			/* as register -disto=file: the file is the next word */
			if (i + 1 >= nb || !word[i + 1]) {
				siril_log_error(_("Missing argument to %s, aborting.\n"), w);
				goto failure;
			}
			g_free(a->disto.filename);
			a->disto.index = DISTO_FILE;
			a->disto.filename = g_canonicalize_filename(word[++i], NULL);	// stored in the file
		}
		else if (!strcmp(w, "-noexpand")) a->framing = FRAMING_CURRENT;	// same as -framing=current
		else if (!strcmp(w, "-expand_force")) a->expand_force = TRUE;
		else if (g_str_has_prefix(w, "-expand_step=")) {
			gchar *end;
			gint64 v = g_ascii_strtoll(w + 13, &end, 10);
			if (end == w + 13 || *end || v < 1 || v > 4096) {
				siril_log_error(_("Wrong parameter values. %s is not valid (1 to 4096), aborting.\n"), w);
				goto failure;
			}
			a->expand_step = (int) v;
		}
		else if (!strcmp(w, "-backup")) a->backup = TRUE;
		else if (!strcmp(w, "-overwrite")) a->overwrite = TRUE;
		else if (g_str_has_prefix(w, "-out=") && w[5]) {
			g_free(a->out);
			a->out = g_strdup(w + 5);
		}
		else {
			siril_log_error(_("Unknown parameter %s, aborting.\n"), w);
			goto failure;
		}
	}
	if (a->overwrite && !a->out) {
		siril_log_error(_("-overwrite only applies to the file given with -out=\n"));
		goto failure;
	}
	if (a->out && lsf_same_file(a->out, a->path)) {
		g_free(a->out);	// the file itself: an append in place
		a->out = NULL;
		a->overwrite = FALSE;
	}
	if (a->out && a->backup) {
		siril_log_error(_("-backup keeps the previous file when it is rewritten; with -out= it is not rewritten\n"));
		goto failure;
	}
	if (a->noreg && a->force_reg) {
		siril_log_error(_("-noreg and -reg cannot be used together\n"));
		goto failure;
	}
	if (a->wcsreg && (a->noreg || a->force_reg)) {
		siril_log_error(_("-wcsreg registers the frames through their WCS and cannot be used with -noreg or -reg\n"));
		goto failure;
	}
	if (a->disto.index != DISTO_UNDEF && a->noreg) {
		siril_log_error(_("-disto undistorts the frames when they are registered and cannot be used with -noreg\n"));
		goto failure;
	}
	if (a->disto.index == DISTO_FILES && !a->wcsreg) {
		siril_log_error(_("-disto=frames undistorts each frame with the SIP terms of its own WCS, which needs -wcsreg\n"));
		goto failure;
	}
	if (a->disto.index == DISTO_FILE) {
		fits dfit = { 0 };
		gchar *msg = NULL;
		gboolean valid = !read_fits_metadata_from_path_first_HDU(a->disto.filename, &dfit)
			&& validate_disto_params(&dfit, a->disto.filename, DISTO_FILE, &msg, NULL);
		clearfits(&dfit);
		if (!valid) {
			siril_log_error(_("LSF: %s cannot be used as a distortion model%s%s\n"), a->disto.filename,
					msg ? ": " : "", msg ? msg : "");
			g_free(msg);
			goto failure;
		}
		g_free(msg);
	}
	if (a->framing == FRAMING_MIN || a->framing == FRAMING_COG) {
		siril_log_error(_("LSF: lsf_append supports -framing=max or -framing=current only. The grid of an LSF "
					"file never shrinks or moves (Section 10.5), so min and cog framing apply to lsf_init\n"));
		goto failure;
	}
	if (lsf_start(a))
		return CMD_GENERIC_ERROR;
	return CMD_OK;
failure:
	lsf_cmd_args_free(a);
	return CMD_ARG_ERROR;
}

/* lsf_export file.lsf.fits output.fit */
int process_lsf_export(int nb) {
	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_EXPORT);
	if (!a)
		return CMD_ALLOC_ERROR;
	a->path = g_strdup(word[1]);
	a->out = g_strdup(word[2]);
	if (lsf_start(a))
		return CMD_GENERIC_ERROR;
	return CMD_OK;
}

/* lsf_util file.lsf.fits [-check] [-repair | -adopt] [-compress[=gzip2|gzip1] | -decompress] [-wcs=file]
 *   [-out=file.lsf.fits [-overwrite]] [-backup] */
int process_lsf_util(int nb) {
	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_UTIL);
	if (!a)
		return CMD_ALLOC_ERROR;
	a->path = g_strdup(word[1]);
	for (int i = 2; i < nb && word[i]; i++) {
		const char *w = word[i];
		int c = parse_compress_option(a, w);
		if (c < 0) goto failure;
		if (c > 0) continue;
		if (!strcmp(w, "-check")) a->check = TRUE;
		else if (!strcmp(w, "-info")) a->info = TRUE;
		else if (!strcmp(w, "-repair")) a->repair = TRUE;
		else if (!strcmp(w, "-adopt")) a->adopt = TRUE;
		else if (!strcmp(w, "-decompress")) {
			if (set_compression_opt(a, LSF_COMP_NONE) < 0) goto failure;
		}
		else if (!strcmp(w, "-backup")) a->backup = TRUE;
		else if (!strcmp(w, "-overwrite")) a->overwrite = TRUE;
		else if (g_str_has_prefix(w, "-out=") && w[5]) {
			g_free(a->out);
			a->out = g_strdup(w + 5);
		}
		else if (g_str_has_prefix(w, "-wcs=") && w[5]) {
			g_free(a->wcs_file);
			a->wcs_file = g_strdup(w + 5);
		}
		else {
			siril_log_error(_("Unknown parameter %s, aborting.\n"), w);
			goto failure;
		}
	}
	if (a->out && lsf_same_file(a->out, a->path)) {
		g_free(a->out);	// the file itself: an update in place
		a->out = NULL;
	}
	if (!a->check && !a->info && !a->repair && !a->adopt && !a->compression_set && !a->wcs_file && !a->out) {
		siril_log_error(_("lsf_util needs an action: -info, -check, -repair, -adopt, -compress, -decompress, -wcs= "
					"or -out=\n"));
		goto failure;
	}
	if (a->info && (a->check || a->repair || a->adopt || a->compression_set || a->wcs_file || a->out)) {
		siril_log_error(_("-info only reads the file and cannot be used with other actions\n"));
		goto failure;
	}
	if (a->repair && a->adopt) {
		siril_log_error(_("-repair and -adopt cannot be used together\n"));
		goto failure;
	}
	if (lsf_start(a))
		return CMD_GENERIC_ERROR;
	return CMD_OK;
failure:
	lsf_cmd_args_free(a);
	return CMD_ARG_ERROR;
}

/* lsf_maps file.lsf.fits prefix [-neff] [-noise] [-sigma] [-snr] [-rejlo] [-rejhi] [-exp] */
int process_lsf_maps(int nb) {
	struct lsf_cmd_args *a = lsf_cmd_args_new(LSF_CMD_MAPS);
	if (!a)
		return CMD_ALLOC_ERROR;
	a->path = g_strdup(word[1]);
	a->out = g_strdup(word[2]);
	for (int i = 3; i < nb && word[i]; i++) {
		const char *w = word[i];
		if (!strcmp(w, "-neff")) a->maps |= LSF_MAP_NEFF;
		else if (!strcmp(w, "-noise")) a->maps |= LSF_MAP_NOISE;
		else if (!strcmp(w, "-sigma")) a->maps |= LSF_MAP_SIGMA;
		else if (!strcmp(w, "-snr")) a->maps |= LSF_MAP_SNR;
		else if (!strcmp(w, "-rejlo")) a->maps |= LSF_MAP_REJLO;
		else if (!strcmp(w, "-rejhi")) a->maps |= LSF_MAP_REJHI;
		else if (!strcmp(w, "-exp")) a->maps |= LSF_MAP_EXP;
		else {
			siril_log_error(_("Unknown parameter %s, aborting.\n"), w);
			lsf_cmd_args_free(a);
			return CMD_ARG_ERROR;
		}
	}
	if (!a->maps)
		a->maps = LSF_MAP_NEFF | LSF_MAP_NOISE | LSF_MAP_SIGMA | LSF_MAP_SNR;
	if (lsf_start(a))
		return CMD_GENERIC_ERROR;
	return CMD_OK;
}
