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

#include <string.h>
#include <math.h>
#include <glib/gstdio.h>

#include "yyjson.h"
#include "core/siril.h"
#include "core/proto.h"
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "core/command.h"
#include "core/command_line_processor.h"
#include "core/processing.h"
#include "core/processing_thread.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "algos/siril_wcs.h"
#include "algos/astrometry_solver.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "lsf/lsf.h"

typedef enum { LSF_CMD_INIT, LSF_CMD_APPEND, LSF_CMD_EXPORT, LSF_CMD_UTIL, LSF_CMD_MAPS } lsf_cmd_mode;

/* maps of lsf_maps (Section 10.8) */
enum { LSF_MAP_NEFF = 1 << 0, LSF_MAP_NOISE = 1 << 1, LSF_MAP_SIGMA = 1 << 2, LSF_MAP_SNR = 1 << 3,
	LSF_MAP_REJLO = 1 << 4, LSF_MAP_REJHI = 1 << 5 };

struct lsf_cmd_args {
	lsf_cmd_mode mode;
	gchar *seqname;
	gchar *path;
	gchar *out;		// export, or the maps prefix
	lsf_settings settings;	// init: requested settings
	gboolean hist_low_set, hist_high_set, hist_margin_set;
	double hist_low, hist_high, hist_margin;
	gboolean noreg, backup, overwrite;
	gboolean rejmaps;		// init and append: -rejmaps (Section 5.4)
	framing_type framing;		// init and append: -framing=, FRAMING_MAX by default (Sections 10.1, 10.5)
	gboolean framing_set;
	gboolean expand_force;
	int expand_step;
	gboolean check, repair, adopt;	// util
	lsf_compression compression;	// -compress=, -nocompress, -decompress (Section 5.3)
	gboolean compression_set;	// append: keeps the compression of the file when not set
	guint maps;		// maps, LSF_MAP_* flags
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

#define LSF_REG_PREFIX "lsf_r_"
static sequence *apply_init_registration(sequence *seq, const struct lsf_cmd_args *a, int *ref,
		int origin[2], Homography *base);
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
	/* registration computed but not applied: lsf_init applies it, by default
	 * onto the union of all the frames (Section 10.1) */
	int ref = sequence_find_refimage(seq);
	sequence *batch = seq;
	int batch_ref = ref;
	if (layer_has_usable_registration(seq, get_registration_layer(seq))) {
		rseq = apply_init_registration(seq, a, &ref, origin, &base);
		if (!rseq || !(src_index = map_to_sources(seq, rseq)))
			goto end;
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
	if (lsf_process_batch(&bargs))
		goto end;
	if (rseq)
		primary.history = g_slist_append(primary.history, g_strdup_printf(
					"LSF gen 1: registration applied, %s framing, origin %d %d", bargs.framing,
					origin[0], origin[1]));
	if (lsf_build_primary(&primary, s, &h))
		goto end;
	uid = g_uuid_string_random();
	if (lsf_write_file(a->path, NULL, s, &h, &primary, uid, 1, rows, 1, FALSE, a->compression))
		goto end;
	siril_log_info(_("LSF file %s created: %u frames, generation 1\n"),
			a->path, rows->len);
	retval = 0;
end:
	g_free(uid);
	clearfits(&primary);
	lsf_history_free(&h);
	g_ptr_array_free(rows, TRUE);
	if (rseq) {
		remove_prefixed_sequence_files(seq, LSF_REG_PREFIX);
		free_sequence(rseq, TRUE);
	}
	free(src_index);
	free_sequence(seq, TRUE);
	return retval;
}

/* registration arguments for Siril's star alignment and apply steps */
static struct registration_args *new_regargs(sequence *seq, const lsf_settings *s) {
	struct registration_args *regargs = calloc(1, sizeof(struct registration_args));
	regargs->seq = seq;
	regargs->reference_image = sequence_find_refimage(seq);
	regargs->layer = seq->nb_layers == 3 ? 1 : 0;
	regargs->filters.filter_included = TRUE;
	regargs->run_in_thread = FALSE;
	regargs->output_scale = 1.f;
	regargs->prefix = strdup(LSF_REG_PREFIX);
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
					"to grow the grid anyway, or -framing=current to crop the batch to the grid\n"),
				f->rx, f->ry, nrx, nry);
		return 1;
	}
	guint64 bytes = (guint64) nrx * nry * f->nchan * (LSF_NB_PL * sizeof(float)
			+ (f->settings.rejmaps ? 3 : 1) * sizeof(WORD));
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

/* registers the batch against HDU 0 of the LSF file (Section 10.4) in two
 * passes: Siril's global star alignment with an external reference computes
 * the transforms, then the footprint of the registered frames decides the
 * expansion of the grid (pad), and the registration is applied onto the
 * grid, expanded or not. Returns the registered sequence, written with the
 * lsf_r_ prefix. The transforms in seq stay relative to the grid before the
 * expansion. */
static sequence *register_against_stack(sequence *seq, const lsf_file *f, const struct lsf_cmd_args *a, int pad[4]) {
	const lsf_settings *s = &f->settings;
	/* the star detection settings stored in the file, over a copy of the
	 * preferences for what they do not hold; the preferences are not changed */
	star_finder_params sf = com.pref.starfinder_conf;
	findstar_from_json(s->findstar_json, &sf);

	/* the registration against the Stack stays in memory and in the
	 * temporary lsf_r_ sequence: the .seq file of the batch and its cached
	 * star lists are left as the user made them */
	struct registration_args *regargs = new_regargs(seq, s);
	regargs->func = register_star_alignment;
	regargs->no_output = TRUE;
	regargs->no_starlist = TRUE;
	regargs->sf_params = &sf;
	regargs->min_pairs = s->min_pairs;
	regargs->max_stars_candidates = MAX_STARS_FITTED;
	regargs->external_ref_path = g_strdup(f->path);
	regargs->use_external_ref = TRUE;

	siril_log_message(_("LSF: registering %d frames against the Stack of %s\n"), seq->selnum, f->path);
	int ret = register_star_alignment(regargs);
	int layer = regargs->layer;
	free_regargs(regargs);
	if (ret || !seq->regparam[layer]) {
		siril_log_error(_("LSF: registration against the Stack failed\n"));
		return NULL;
	}

	/* footprint of the registered frames on the grid */
	double box[4];
	int first = -1;
	lsf_footprint_init(box);
	for (int i = 0; i < seq->number; i++) {
		const Homography *H = &seq->regparam[layer][i].H;
		if (!seq->imgparam[i].incl || guess_transform_from_H(*H) == NULL_TRANSFORMATION)
			continue;
		double h[9] = { H->h00, H->h01, H->h02, H->h10, H->h11, H->h12, H->h20, H->h21, H->h22 };
		lsf_footprint_add(h, seq->is_variable ? seq->imgparam[i].rx : seq->rx,
				seq->is_variable ? seq->imgparam[i].ry : seq->ry, box);
		if (first < 0)
			first = i;
	}
	if (first < 0) {
		siril_log_error(_("LSF: no frame could be registered against the Stack\n"));
		return NULL;
	}
	if (plan_expansion(f, a, box, pad))
		return NULL;

	/* second pass. Apply checks the reference frame of the sequence, which
	 * plays no part with an external reference: use a registered one. */
	regargs = new_regargs(seq, s);
	regargs->func = register_apply_reg;
	regargs->no_output = FALSE;
	if (any_pad(pad))
		regargs->canvas =(framing_roi) { -pad[LSF_PAD_LEFT], -pad[LSF_PAD_LAST],
			f->rx + pad[LSF_PAD_LEFT] + pad[LSF_PAD_RIGHT], f->ry + pad[LSF_PAD_FIRST] + pad[LSF_PAD_LAST] };
	int seq_ref = seq->reference_image;
	seq->reference_image = first;
	regargs->reference_image = first;
	ret = register_apply_reg(regargs);
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
static sequence *apply_init_registration(sequence *seq, const struct lsf_cmd_args *a, int *ref,
		int origin[2], Homography *base) {
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
	Homography Href = seq->regparam[layer][*ref].H;
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
			cvTransfH(&Hi, &Href, &Hr);	// Href^-1 * Hi
			double h[9] = { Hr.h00, Hr.h01, Hr.h02, Hr.h10, Hr.h11, Hr.h12, Hr.h20, Hr.h21, Hr.h22 };
			lsf_footprint_add(h, seq->is_variable ? seq->imgparam[i].rx : seq->rx,
					seq->is_variable ? seq->imgparam[i].ry : seq->ry, box);
		}
		lsf_grid_overhang(box, rx, ry, pad);
	}
	struct registration_args *regargs = new_regargs(seq, &a->settings);
	regargs->func = register_apply_reg;
	regargs->layer = layer;
	regargs->no_output = FALSE;
	regargs->framing = a->framing == FRAMING_MAX ? FRAMING_CURRENT : a->framing;
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
	free_regargs(regargs);
	*base = Ht;
	cvInvertH(base);
	cvMultH(*base, Href, &R);
	origin[0] = (int) lround(R.h01 * (ry - 1) + R.h02);
	origin[1] = nry - 1 - (int) lround(R.h11 * (ry - 1) + R.h12);
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
	GPtrArray *rows = g_ptr_array_new_with_free_func(lsf_ledger_row_free);

	if (lsf_open_and_validate(a->path, &f, FALSE)) {
		lsf_file_free(&f);
		g_ptr_array_free(rows, TRUE);
		return 1;
	}
	/* this writer merges onto a grid at the frames' own scale only: frames
	 * merged onto a drizzled grid (Section 6) would land on the wrong pixels */
	if (f.settings.oversample != 1) {
		siril_log_error(_("LSF: %s has a grid drizzled at scale %d, which this version cannot append to. "
					"Only grid.oversample 1 is supported\n"), a->path, f.settings.oversample);
		goto end;
	}
	double hist_before[3] = { f.settings.hist_low, f.settings.hist_high, f.settings.hist_margin };
	apply_hist_overrides(a, &f.settings);
	if (hist_before[0] != f.settings.hist_low || hist_before[1] != f.settings.hist_high
			|| hist_before[2] != f.settings.hist_margin)
		siril_log_message(_("LSF: history test changed from low %g high %g margin %g to low %g high %g margin %g. "
					"It applies from this append on and is saved in the file for the following ones; "
					"the frames already in the stack are not affected\n"),
				hist_before[0], hist_before[1], hist_before[2],
				f.settings.hist_low, f.settings.hist_high, f.settings.hist_margin);
	if (a->rejmaps && f.settings.rejmaps)
		siril_log_message(_("LSF: %s already has rejection maps, counted since generation %d\n"),
				a->path, f.hist.rej_gen);
	f.settings.rejmaps |= a->rejmaps;	// lsf_process_batch() creates them
	seq = load_batch(a->seqname);
	if (!seq)
		goto end;
	if (seq->nb_layers != f.nchan) {
		siril_log_error(_("LSF: the file has %d channel(s), the batch has %d\n"), f.nchan, seq->nb_layers);
		goto end;
	}

	sequence *batch = seq;
	int gen = f.gen + 1, batchnum = f.nb_batches + 1;
	int pad[4] = { 0 };
	if (!a->noreg) {
		rseq = register_against_stack(seq, &f, a, pad);
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
	if (lsf_write_file(a->path, &f, &f.settings, &f.hist, &f.primary, f.uid, gen, rows, batchnum, a->backup, comp))
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

static void print_summary(const lsf_file *f) {
	const lsf_settings *s = &f->settings;
	const fkeywords *k = &f->primary.keywords;
	siril_log_message(_("LSF file %s: version %s, id %s, generation %d\n"), f->path, LSF_VERSION, f->uid, f->gen);
	siril_log_message(_("  %ld frames in %d batch(es), STACKCNT %u, live time %.1f s\n"), f->ledger_rows,
			f->nb_batches, k->stackcnt, k->livetime > 0.0 ? k->livetime : 0.0);
	siril_log_message(_("  grid %dx%d, %d channel(s), sensor %s%s%s\n"), f->rx, f->ry, f->nchan,
			s->is_cfa ? "CFA " : "mono", s->is_cfa ? s->bayer : "", s->is_cfa ? " (debayered)" : "");
	siril_log_message(_("  rejection %s %g %g, history test low %g high %g margin %g\n"),
			lsf_rejection_name(s->rej), s->sig[0], s->sig[1], s->hist_low, s->hist_high, s->hist_margin);
	siril_log_message(_("  normalization %s%s%s, output normalization %s\n"), lsf_norm_name(s->norm),
			s->lite_norm ? " (fast)" : "", s->rgb_equal ? " (RGB equalized)" : "", s->outnorm ? "on" : "off");
	switch (s->weighting) {
	case NOISE_WEIGHT:
		siril_log_message(_("  weighting noise, frozen norm %g %g %g\n"), s->noise_norm[0],
				f->nchan == 3 ? s->noise_norm[1] : 0.0, f->nchan == 3 ? s->noise_norm[2] : 0.0);
		break;
	case WFWHM_WEIGHT:
		siril_log_message(_("  weighting wfwhm, frozen min %g max %g norm %g\n"), s->fwhm_min, s->fwhm_max, s->fwhm_norm);
		break;
	case NBSTARS_WEIGHT:
		siril_log_message(_("  weighting nbstars, frozen min %g max %g norm %g\n"), s->stars_min, s->stars_max, s->stars_norm);
		break;
	default:
		siril_log_message(_("  weighting none\n"));
	}
	double raw = (double) f->rx * f->ry * f->nchan;
	siril_log_message(_("  LSF_ACC %s, %.0f MB stored (%.0f%% of raw), LSF_COUNT %s, %.0f MB stored (%.0f%% of raw)\n"),
			f->acc_codec, (double) f->acc_bytes / BYTES_IN_A_MB, 100.0 * f->acc_bytes / (raw * LSF_NB_PL * sizeof(float)),
			f->count_codec, (double) f->count_bytes / BYTES_IN_A_MB, 100.0 * f->count_bytes / (raw * sizeof(WORD)));
	if (f->hist.rej[0][0])
		siril_log_message(_("  LSF_REJ %s, %.0f MB stored (%.0f%% of raw), rejections counted since generation %d%s\n"),
				f->rej_codec, (double) f->rej_bytes / BYTES_IN_A_MB, 100.0 * f->rej_bytes / (raw * 2 * sizeof(WORD)),
				f->hist.rej_gen, f->hist.rej_saturated ? ", some counts saturated" : "");
	else siril_log_message(_("  no rejection maps (LSF_REJ)\n"));
}

/* lsf_util: Section 8 validation and a summary (-check), Repair or Adopt
 * when only check 8 failed, and a change of the compression of the
 * extensions (Section 5.3). The actions combine into one rewrite, which
 * increments the generation. */
static int lsf_util_run(struct lsf_cmd_args *a) {
	lsf_file f;
	int retval = 1;
	GPtrArray *none = NULL;
	if (lsf_open_and_validate(a->path, &f, TRUE)) {
		lsf_file_free(&f);
		return 1;
	}
	print_summary(&f);
	if (!f.stack_modified)
		siril_log_info(_("LSF: %s passed all the checks of Section 8\n"), a->path);
	gboolean fix_stack = (a->repair || a->adopt) && f.stack_modified;
	if ((a->repair || a->adopt) && !f.stack_modified)
		siril_log_message(_("LSF: the Stack matches the accumulator, there is nothing to repair or adopt\n"));
	gboolean recompress = a->compression_set && a->compression != f.compression;
	if (a->compression_set && !recompress) {
		if (f.compression == LSF_COMP_NONE)
			siril_log_message(_("LSF: the extensions of %s are already uncompressed\n"), a->path);
		else siril_log_message(_("LSF: the extensions of %s are already compressed with %s\n"), a->path,
				lsf_compression_name(f.compression));
	}
	/* any rewrite rebuilds the Stack from MEAN, which would silently repair it */
	if (f.stack_modified && !fix_stack) {
		siril_log_error(_("LSF: %s fails check 8, its Stack was modified outside Siril. lsf_append will "
					"refuse it, and lsf_util will not rewrite it, until lsf_util is run with -repair (rebuild "
					"the Stack from the accumulator) or -adopt (keep the edit)\n"), a->path);
		goto end;
	}
	/* -out writes the result to a new file and leaves the source untouched,
	 * even when no action changes anything */
	const char *target = a->out ? a->out : a->path;
	if (!fix_stack && !recompress && !a->out) {
		retval = 0;
		goto end;
	}
	if (a->out && g_file_test(a->out, G_FILE_TEST_EXISTS) && !a->overwrite) {
		siril_log_error(_("LSF: %s already exists, use -overwrite to replace it\n"), a->out);
		goto end;
	}
	gchar *card;
	int gen = f.gen + 1;
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
	if (lsf_build_primary(&f.primary, &f.settings, &f.hist))
		goto end;
	none = g_ptr_array_new();
	if (lsf_write_file(target, &f, &f.settings, &f.hist, &f.primary, f.uid, gen, none, f.nb_batches, a->backup, comp))
		goto end;
	if (fix_stack && a->adopt)
		siril_log_info(_("LSF: edited Stack of %s adopted, generation %d\n"), target, gen);
	else if (fix_stack)
		siril_log_info(_("LSF: Stack of %s repaired from the accumulator, generation %d\n"), target, gen);
	if (recompress && comp == LSF_COMP_NONE)
		siril_log_info(_("LSF: extensions of %s decompressed, generation %d\n"), target, gen);
	else if (recompress)
		siril_log_info(_("LSF: extensions of %s compressed with %s, generation %d\n"), target,
				lsf_compression_name(comp), gen);
	if (a->out)
		siril_log_info(_("LSF: %s written from %s, generation %d, %s left unchanged\n"), a->out, a->path,
				gen, a->path);
	retval = 0;
end:
	if (none)
		g_ptr_array_free(none, TRUE);
	lsf_file_free(&f);
	return retval;
}

/* lsf_maps: images derived from the accumulator (Sections 10.7, 10.8). The
 * noise and sigma maps are in the units of the Stack. */
static int lsf_maps_run(struct lsf_cmd_args *a) {
	static const struct { guint flag; int rej; const char *suffix, *what; } kinds[] = {
		{ LSF_MAP_NEFF, -1, "neff", "effective number of frames" },
		{ LSF_MAP_NOISE, -1, "noise", "standard deviation of the mean" },
		{ LSF_MAP_SIGMA, -1, "sigma", "standard deviation of one frame" },
		{ LSF_MAP_SNR, -1, "snr", "mean over its standard deviation" },
		{ LSF_MAP_REJLO, 0, "rejlo", "samples rejected low" },
		{ LSF_MAP_REJHI, 1, "rejhi", "samples rejected high" },
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
		for (int c = 0; rej < 0 && c < f.nchan; c++) {
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
		else fit->history = g_slist_append(fit->history, g_strdup_printf("LSF %s map (%s)",
					kinds[m].suffix, kinds[m].what));
		fit->history = g_slist_append(fit->history, g_strdup_printf("from LSF file %s, generation %d",
					f.uid, f.gen));
		gchar *name = g_strdup_printf("%s_%s", a->out, kinds[m].suffix);
		if (savefits(name, fit)) {
			siril_log_error(_("LSF: could not save the %s map\n"), kinds[m].suffix);
			retval = 1;
		} else siril_log_message(_("LSF: %s map (%s) saved as %s\n"), kinds[m].suffix, kinds[m].what, name);
		g_free(name);
		clearfits(fit);
		free(fit);
	}
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
	case LSF_CMD_UTIL:
		retval = lsf_util_run(a);
		break;
	case LSF_CMD_MAPS:
		retval = lsf_maps_run(a);
		break;
	default:
		retval = lsf_export_master(a->path, a->out);
		if (!retval) {
			/* the name savefits gave the file */
			gchar *name = set_right_extension(a->out);
			siril_log_message(_("LSF: master of %s exported to %s\n"), a->path, name ? name : a->out);
			g_free(name);
		}
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
	/* thresholds in sigma units: 0 would reject every sample on that side of the mean */
	if (end == value || *end || !(v > 0.0) || !isfinite(v)) {
		siril_log_error(_("Wrong parameter values. %s is not valid, aborting.\n"), word_i);
		*set = -1;
		return TRUE;
	}
	*out = v;
	*set = TRUE;
	return TRUE;
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
	if (!strcmp(w, "-nocompress"))
		return set_compression_opt(a, LSF_COMP_NONE);
	return parse_compress_option(a, w);
}

/* lsf_init sequencename file.lsf.fits [rej [type] sigma_low sigma_high] [-nonorm, -norm=] [-fastnorm]
 *   [-rgb_equal] [-weight=] [-output_norm] [-framing=max|current|min|cog] [-hist_low=] [-hist_high=] [-hist_margin=]
 *   [-rejmaps] [-compress[=gzip2|gzip1]] [-overwrite] */
int process_lsf_init(int nb) {
	struct lsf_cmd_args *a = calloc(1, sizeof(struct lsf_cmd_args));
	a->mode = LSF_CMD_INIT;
	a->framing = FRAMING_MAX;
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

/* lsf_append file.lsf.fits sequencename [-noreg] [-framing=max|current] [-noexpand] [-expand_step=] [-expand_force] [-backup]
 *   [-hist_low=] [-hist_high=] [-hist_margin=] [-rejmaps] [-compress[=gzip2|gzip1] | -nocompress] */
int process_lsf_append(int nb) {
	struct lsf_cmd_args *a = calloc(1, sizeof(struct lsf_cmd_args));
	a->mode = LSF_CMD_APPEND;
	a->path = g_strdup(word[1]);
	a->seqname = g_strdup(word[2]);
	a->expand_step = 64;
	a->framing = FRAMING_MAX;
	for (int i = 3; i < nb && word[i]; i++) {
		const char *w = word[i];
		int c = parse_common_option(a, w);
		if (c < 0) goto failure;
		if (c > 0) continue;
		if (!strcmp(w, "-noreg")) a->noreg = TRUE;
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
		else {
			siril_log_error(_("Unknown parameter %s, aborting.\n"), w);
			goto failure;
		}
	}
	if (a->framing == FRAMING_MIN || a->framing == FRAMING_COG) {
		siril_log_error(_("LSF: lsf_append supports -framing=max or -framing=current only. The grid of an LSF "
					"file never shrinks or moves (Section 10.5), so min and cog framing apply to lsf_init\n"));
		goto failure;
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

/* lsf_util file.lsf.fits [-check] [-repair | -adopt] [-compress[=gzip2|gzip1] | -decompress] [-out=file.lsf.fits [-overwrite]] [-backup] */
int process_lsf_util(int nb) {
	struct lsf_cmd_args *a = calloc(1, sizeof(struct lsf_cmd_args));
	a->mode = LSF_CMD_UTIL;
	a->path = g_strdup(word[1]);
	for (int i = 2; i < nb && word[i]; i++) {
		const char *w = word[i];
		int c = parse_compress_option(a, w);
		if (c < 0) goto failure;
		if (c > 0) continue;
		if (!strcmp(w, "-check")) a->check = TRUE;
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
		else {
			siril_log_error(_("Unknown parameter %s, aborting.\n"), w);
			goto failure;
		}
	}
	if (a->out && !strcmp(a->out, a->path)) {
		g_free(a->out);	// the file itself: an update in place
		a->out = NULL;
	}
	if (!a->check && !a->repair && !a->adopt && !a->compression_set && !a->out) {
		siril_log_error(_("lsf_util needs an action: -check, -repair, -adopt, -compress, -decompress or -out=\n"));
		goto failure;
	}
	if (a->repair && a->adopt) {
		siril_log_error(_("-repair and -adopt cannot be used together\n"));
		goto failure;
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

/* lsf_maps file.lsf.fits prefix [-neff] [-noise] [-sigma] [-snr] [-rejlo] [-rejhi] */
int process_lsf_maps(int nb) {
	struct lsf_cmd_args *a = calloc(1, sizeof(struct lsf_cmd_args));
	a->mode = LSF_CMD_MAPS;
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
		else {
			siril_log_error(_("Unknown parameter %s, aborting.\n"), w);
			free_cmd_args(a);
			return CMD_ARG_ERROR;
		}
	}
	if (!a->maps)
		a->maps = LSF_MAP_NEFF | LSF_MAP_NOISE | LSF_MAP_SIGMA | LSF_MAP_SNR;
	if (!start_in_new_thread(lsf_worker, a)) {
		free_cmd_args(a);
		return CMD_GENERIC_ERROR;
	}
	return CMD_OK;
}
