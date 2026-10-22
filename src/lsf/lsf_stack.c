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

/* LSF batch engine (FITS profile v0.4, Sections 10.2, 10.4, 10.7).
 *
 * The batch, a sequence already on the LSF grid, is read in blocks of rows
 * like stack_mean_with_rejection() does. Each pixel goes through Siril's own
 * normalization formulas (with coefficients computed against the frozen
 * baseline of the file), Siril's batch-local rejection
 * (apply_rejection_float), then the LSF history test. The accepted samples
 * are merged into the accumulator with the weighted pairwise update.
 * Siril's stacking code is not modified. */

#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <gsl/gsl_cdf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/siril.h"
#include "core/proto.h"
#include "core/processing.h"
#include "core/siril_log.h"
#include "core/OS_utils.h"
#include "core/gui_iface.h"
#include "algos/statistics.h"
#include "io/sequence.h"
#include "io/image_format_fits.h"
#include "opencv/opencv.h"
#include "registration/registration.h"
#include "lsf/lsf.h"

/* from median_and_mean.c */
int stack_open_all_files(struct stacking_args *args, int *bitpix, int *naxis, long *naxes,
		GList **list_date, fits *fit);

void lsf_ledger_row_free(gpointer p) {
	lsf_ledger_row *r = p;
	if (!r) return;
	g_free(r->frame_id);
	g_free(r->filename);
	g_free(r->date_obs);
	g_free(r->filter);
	g_free(r->orighdr);
	g_free(r->align_json);
	g_free(r);
}

/* per-frame data of the batch */
typedef struct {
	int nb;			// frames kept
	int *idx;		// sequence index of each kept frame
	double *pscale[3], *poffset[3], *pmul[3];
	double *weight[3];
	double *sframe[3];	// frame noise in accumulator units (sigma_frame)
} batch_frames;

static void free_batch_frames(batch_frames *b) {
	free(b->idx);
	for (int c = 0; c < 3; c++) {
		free(b->pscale[c]); free(b->poffset[c]); free(b->pmul[c]);
		free(b->weight[c]); free(b->sframe[c]);
	}
}

static double est_loc(const imstats *st, gboolean lite) { return lite ? st->median : st->location; }
static double est_scale(const imstats *st, gboolean lite) { return lite ? 1.5 * st->mad : st->scale; }

static gchar *format_date(GDateTime *dt) {
	if (!dt) return NULL;
	gchar *base = g_date_time_format(dt, "%Y-%m-%dT%H:%M:%S");
	gchar *s = g_strdup_printf("%s.%03d", base, g_date_time_get_microsecond(dt) / 1000);
	g_free(base);
	return s;
}

static gboolean dbl_known(double v) { return v > DEFAULT_DOUBLE_VALUE; }

static gboolean differ_dbl(double a, double b) {
	return dbl_known(a) && dbl_known(b) && fabs(a - b) > 1e-6 * MAX(fabs(a), fabs(b));
}

static gboolean differ_str(const char *a, const char *b) {
	return a[0] != '\0' && b[0] != '\0' && strcmp(a, b);
}

/* Section 7, keywords set at initialization then checked on append */
static const char *check_frame_keywords(const fits *primary, const fits *frame) {
	const fkeywords *p = &primary->keywords, *f = &frame->keywords;
	if (p->binning_x && f->binning_x && (p->binning_x != f->binning_x || p->binning_y != f->binning_y))
		return "XBINNING/YBINNING";
	if (differ_dbl(p->pixel_size_x, f->pixel_size_x) || differ_dbl(p->pixel_size_y, f->pixel_size_y))
		return "XPIXSZ/YPIXSZ";
	if (differ_dbl(p->focal_length, f->focal_length))
		return "FOCALLEN";
	if (differ_str(p->instrume, f->instrume))
		return "INSTRUME";
	if (differ_str(p->filter, f->filter))
		return "FILTER";
	return NULL;
}

/* Section 7, keywords copied from the first frame and removed when a later
 * frame differs. Siril does not write a keyword holding its default value. */
static void reconcile_first_frame_keywords(fits *primary, const fits *frame) {
	fkeywords *p = &primary->keywords;
	const fkeywords *f = &frame->keywords;
	if (differ_dbl(p->exposure, f->exposure)) p->exposure = DEFAULT_DOUBLE_VALUE;
	if (p->key_gain != (int) DEFAULT_UINT_VALUE && f->key_gain != (int) DEFAULT_UINT_VALUE && p->key_gain != f->key_gain)
		p->key_gain = DEFAULT_UINT_VALUE;
	if (p->key_offset != (int) DEFAULT_UINT_VALUE && f->key_offset != (int) DEFAULT_UINT_VALUE && p->key_offset != f->key_offset)
		p->key_offset = DEFAULT_UINT_VALUE;
	if (differ_dbl(p->ccd_temp, f->ccd_temp)) p->ccd_temp = DEFAULT_DOUBLE_VALUE;
	if (differ_dbl(p->cvf, f->cvf)) p->cvf = DEFAULT_DOUBLE_VALUE;
	if (differ_str(p->telescop, f->telescop)) p->telescop[0] = '\0';
	if (differ_str(p->object, f->object)) p->object[0] = '\0';
}

/* Siril keeps a frame's header as text, one card per line. FITS string columns
 * hold printable ASCII only, so store it the FITS way: 80-character cards one
 * after another, no line breaks (what fits_hdr2str produces). */
static gchar *header_as_cards(const char *header) {
	if (!header)
		return g_strdup("");
	GString *out = g_string_sized_new(strlen(header) + 80);
	gchar **lines = g_strsplit(header, "\n", -1);
	for (int i = 0; lines[i]; i++) {
		char card[81];
		int n = 0;
		for (const char *c = lines[i]; *c && n < 80; c++)
			if (*c >= 32 && *c <= 126)
				card[n++] = *c;
		if (n == 0)
			continue;
		memset(card + n, ' ', 80 - n);
		g_string_append_len(out, card, 80);
	}
	g_strfreev(lines);
	return g_string_free(out, FALSE);
}

static void fill_ledger_row(lsf_ledger_row *r, struct lsf_batch_args *args, const batch_frames *b,
		int k, const fits *meta, imstats **stats) {
	sequence *seq = args->seq;
	int idx = b->idx[k];
	int nchan = seq->nb_layers;
	int reglayer = nchan == 3 ? 1 : 0;
	sequence *src = args->src_seq ? args->src_seq : seq;
	int src_idx = args->src_index ? args->src_index[idx] : idx;
	char name[256];

	r->frame_id = g_strdup_printf("%.20s_%05d", args->seq_name, src->imgparam[src_idx].filenum);
	if (seq_get_image_filename(src, src_idx, name))
		r->filename = g_strdup(name);
	r->date_obs = format_date(meta->keywords.date_obs);
	r->exptime = dbl_known(meta->keywords.exposure) ? meta->keywords.exposure : NAN;
	r->gain = meta->keywords.key_gain != (int) DEFAULT_UINT_VALUE ? meta->keywords.key_gain : INT_MIN;
	r->offset = meta->keywords.key_offset != (int) DEFAULT_UINT_VALUE ? meta->keywords.key_offset : INT_MIN;
	r->xbin = meta->keywords.binning_x ? (int) meta->keywords.binning_x : INT_MIN;
	r->ybin = meta->keywords.binning_y ? (int) meta->keywords.binning_y : INT_MIN;
	r->ccd_temp = dbl_known(meta->keywords.ccd_temp) ? (float) meta->keywords.ccd_temp : NAN;
	r->filter = g_strdup(meta->keywords.filter);
	r->orighdr = header_as_cards(meta->header);

	/* registration: the transform onto the grid is the one computed for the
	 * source frame when the batch was registered against the Stack, or, when
	 * lsf_init applied the registration of the sequence, that one mapped onto
	 * the grid */
	int src_layer = args->xfrm_base ? get_registration_layer(src) : reglayer;
	regdata *src_reg = (src_layer >= 0 && layer_has_registration(src, src_layer)) ? &src->regparam[src_layer][src_idx] : NULL;
	regdata *reg = (layer_has_registration(seq, reglayer)) ? &seq->regparam[reglayer][idx] : NULL;
	Homography H;
	if (args->src_seq && src_reg) {
		H = src_reg->H;
		if (args->xfrm_base)
			cvMultH(*args->xfrm_base, src_reg->H, &H);
	} else cvGetEye(&H);
	double h[9] = { H.h00, H.h01, H.h02, H.h10, H.h11, H.h12, H.h20, H.h21, H.h22 };
	memcpy(r->xfrm, h, sizeof h);
	static const char *transf[] = { "SHIFT", "SIMILAR", "AFFINE", "HOMOGR" };
	transformation_type t = args->settings->transform;
	if (args->xfrm_base && src_reg) {
		t = guess_transform_from_H(src_reg->H);	// unknown otherwise: the user registered the sequence
		if (t < SHIFT_TRANSFORMATION) t = SHIFT_TRANSFORMATION;
	}
	g_strlcpy(r->reg_type, args->src_seq ? transf[t < 0 ? 3 : t] : "NONE", sizeof r->reg_type);
	if (!args->src_seq)
		r->align_json = g_strdup("{\"reference\":\"none\"}");
	else if (args->xfrm_base)
		r->align_json = g_strdup_printf("{\"reference\":\"sequence\",\"framing\":\"%s\",\"transform\":\"%s\",\"interpolation\":%d}",
				args->framing ? args->framing : "current", transf[t], args->settings->interpolation);
	else r->align_json = g_strdup_printf("{\"reference\":\"stack\",\"framing\":\"%s\",\"transform\":\"%s\",\"interpolation\":%d,\"min_pairs\":%d}",
				args->framing ? args->framing : "current", transf[t < 0 ? 3 : t], args->settings->interpolation,
				args->settings->min_pairs);

	regdata *q = reg ? reg : src_reg;
	r->fwhm = q ? q->fwhm : NAN;
	r->wfwhm = q ? q->weighted_fwhm : NAN;
	r->round = q ? q->roundness : NAN;
	r->regq = q ? (float) q->quality : NAN;
	r->regbg = q ? q->background_lvl : NAN;
	r->nstars = q ? q->number_of_stars : INT_MIN;

	for (int c = 0; c < nchan; c++) {
		r->norm_sc[c] = b->pscale[c][k] * b->pmul[c][k];
		r->norm_off[c] = b->poffset[c][k];
		g_strlcpy(r->norm_mth[c], lsf_norm_name(args->settings->norm), sizeof r->norm_mth[c]);
		r->weight[c] = b->weight[c][k];
		imstats *st = stats[c];
		r->bkg_med[c] = st ? (float) st->median : NAN;
		r->noise[c] = st ? (float) st->bgnoise : NAN;
		r->noise_n[c] = (float) b->sframe[c][k];
		r->st_mean[c] = st ? (float) st->mean : NAN;
		r->st_sig[c] = st ? (float) st->sigma : NAN;
		r->st_mad[c] = st ? (float) st->mad : NAN;
		r->st_sqbw[c] = st ? (float) st->sqrtbwmv : NAN;
		r->st_loc[c] = st ? (float) st->location : NAN;
		r->st_scl[c] = st ? (float) st->scale : NAN;
	}
}

/* per-thread buffers, same layout needs as stack_mean_with_rejection() */
typedef struct {
	struct _data_block db;
	float **pix;		// [frame] block of rows
	WORD *wbuf;		// ushort read buffer
	float *xbuf;		// accepted samples
	double *wbufd;		// their weights
	gboolean *sel;		// frame survived the batch-local rejection
	void *mem;
} lsf_thread_data;

static int alloc_thread_data(lsf_thread_data *t, int nb, size_t block_pix, gboolean ushort_input, rejection rej) {
	memset(t, 0, sizeof(lsf_thread_data));
	t->pix = malloc(nb * sizeof(float *));
	t->mem = malloc((size_t) nb * block_pix * sizeof(float));
	t->db.stack = malloc(nb * sizeof(float));
	t->db.o_stack = malloc(nb * sizeof(float));
	t->db.w_stack = malloc(nb * sizeof(float));
	t->db.rejected = malloc(nb * sizeof(int));
	t->xbuf = malloc(nb * sizeof(float));
	t->wbufd = malloc(nb * sizeof(double));
	t->sel = malloc(nb * sizeof(gboolean));
	if (ushort_input)
		t->wbuf = malloc(block_pix * sizeof(WORD));
	if (!t->pix || !t->mem || !t->db.stack || !t->db.o_stack || !t->db.w_stack || !t->db.rejected
			|| !t->xbuf || !t->wbufd || !t->sel || (ushort_input && !t->wbuf))
		return 1;
	for (int f = 0; f < nb; f++)
		t->pix[f] = (float *) t->mem + (size_t) f * block_pix;
	if (rej == LINEARFIT) {
		t->db.xf = malloc(2 * nb * sizeof(float));
		if (!t->db.xf) return 1;
		t->db.yf = t->db.xf + nb;
		/* precomputation, as in stack_mean_or_median() */
		t->db.m_x = (nb - 1) * 0.5f;
		t->db.m_dx2 = 0.f;
		for (int j = 0; j < nb; ++j) {
			const float dx = j - t->db.m_x;
			t->db.xf[j] = 1.f / (j + 1);
			t->db.m_dx2 += (dx * dx - t->db.m_dx2) * t->db.xf[j];
		}
		t->db.m_dx2 = 1.f / t->db.m_dx2;
	}
	return 0;
}

static void free_thread_data(lsf_thread_data *t) {
	free(t->pix); free(t->mem); free(t->wbuf); free(t->xbuf); free(t->wbufd); free(t->sel);
	free(t->db.stack); free(t->db.o_stack); free(t->db.w_stack); free(t->db.rejected);
	free(t->db.xf);
}

static double normalize_sample(normalization norm, float v, double pscale, double poffset, double pmul) {
	switch (norm) {
	case ADDITIVE:
	case ADDITIVE_SCALING:
		return v != 0.f ? v * pscale - poffset : 0.0;
	case MULTIPLICATIVE:
	case MULTIPLICATIVE_SCALING:
		return v * pscale * pmul;
	default:
		return v;
	}
}

int lsf_process_batch(struct lsf_batch_args *args) {
	sequence *seq = args->seq;
	lsf_settings *s = args->settings;
	lsf_history *h = args->hist;
	int nchan = seq->nb_layers;
	int reglayer = nchan == 3 ? 1 : 0;
	int retval = 1;
	batch_frames b = { 0 };
	fits *metas = NULL;
	imstats ***stats = NULL;	// [frame][chan]
	int nb_incl = 0;
	struct stacking_args sargs = { 0 };
	GList *list_date = NULL;
	lsf_thread_data *tdata = NULL;
	struct _image_block *blocks = NULL;
	int nb_threads = 1;
	gboolean opened = FALSE;
	fits refmeta = { 0 };
	struct stacking_args rargs = { 0 };	// rejection arguments for apply_rejection_float()

	if (seq->is_variable) {
		siril_log_error(_("LSF: the frames of the batch must all have the size of the grid\n"));
		return 1;
	}
	if (nchan != 1 && nchan != 3) {
		siril_log_error(_("LSF: unsupported number of channels (%d)\n"), nchan);
		return 1;
	}
	if (h->plane[0][0]) {
		if (seq->rx != h->rx || seq->ry != h->ry || nchan != h->nchan) {
			siril_log_error(_("LSF: the batch is %dx%dx%d, the LSF grid is %dx%dx%d. Frames that are "
						"not on the grid must be registered by lsf_append (no -noreg).\n"),
					seq->rx, seq->ry, nchan, h->rx, h->ry, h->nchan);
			return 1;
		}
	} else if (lsf_history_alloc(h, seq->rx, seq->ry, nchan)) {
		return 1;
	}

	/* 1. included frames, their metadata and the Section 7 checks */
	int *cand = malloc(seq->number * sizeof(int));
	metas = calloc(seq->number, sizeof(fits));
	stats = calloc(seq->number, sizeof(imstats **));
	if (!cand || !metas || !stats) {
		PRINT_ALLOC_ERR;
		free(cand);
		goto end;
	}
	for (int i = 0; i < seq->number; i++) {
		if (!seq->imgparam[i].incl)
			continue;
		if (seq_read_frame_metadata(seq, i, &metas[i])) {
			free(cand);
			goto end;
		}
		if (!args->init) {
			const char *bad = check_frame_keywords(args->primary, &metas[i]);
			if (bad) {
				siril_log_warning(_("LSF: frame %d not merged, %s differs from the LSF file\n"),
						seq->imgparam[i].filenum, bad);
				continue;
			}
		}
		cand[nb_incl++] = i;
	}
	if (args->init) {
		gboolean found = FALSE;
		for (int k = 0; k < nb_incl; k++)
			found |= cand[k] == args->ref_image;
		if (!found && nb_incl > 0) {
			/* like stack (sequence_filtering.c): the first selected frame is
			 * the temporary reference for this run */
			siril_log_warning(_("LSF: the reference image is not in the selected set of images, "
					"using image %d as temporary reference image\n"), cand[0] + 1);
			args->ref_image = cand[0];
		}
	}
	if (nb_incl == 0) {
		siril_log_error(_("LSF: no frame to merge in this batch\n"));
		free(cand);
		goto end;
	}

	/* 2. statistics, Siril's estimators (cached in the sequence) */
	gui_iface.set_progress(PROGRESS_RESET, _("LSF: computing frame statistics"));
	for (int k = 0; k < nb_incl; k++) {
		int i = cand[k];
		stats[i] = calloc(3, sizeof(imstats *));
		if (compute_all_channels_statistics_seqimage(seq, i, NULL, s->lite_norm ? STATS_LITENORM : STATS_NORM,
					MULTI_THREADED, -1, stats[i])) {
			siril_log_error(_("LSF: could not compute the statistics of frame %d\n"), seq->imgparam[i].filenum);
			free(cand);
			goto end;
		}
		for (int c = 0; c < nchan; c++)
			if (!stats[i][c]) {
				free(cand);
				goto end;
			}
		gui_iface.set_progress((double) (k + 1) / nb_incl, NULL);
		if (!processing_should_continue()) {
			free(cand);
			goto end;
		}
	}
	writeseqfile(seq);

	/* 3. normalization against the frozen baseline */
	if (args->init && s->norm != NO_NORM) {
		imstats **rs = stats[args->ref_image];
		for (int c = 0; c < nchan; c++) {
			s->base_loc[c] = est_loc(rs[c], s->lite_norm);
			s->base_scale[c] = est_scale(rs[c], s->lite_norm);
		}
	}
	int nb = nb_incl;
	b.idx = malloc(nb * sizeof(int));
	for (int c = 0; c < nchan; c++) {
		b.pscale[c] = malloc(nb * sizeof(double));
		b.poffset[c] = malloc(nb * sizeof(double));
		b.pmul[c] = malloc(nb * sizeof(double));
		b.weight[c] = malloc(nb * sizeof(double));
		b.sframe[c] = malloc(nb * sizeof(double));
	}
	double *noise_raw = malloc(nchan * nb * sizeof(double));
	float *wfwhm = malloc(nb * sizeof(float));
	int *nstars = malloc(nb * sizeof(int));
	gboolean has_reg = layer_has_registration(seq, reglayer);
	for (int k = 0; k < nb; k++) {
		int i = cand[k];
		for (int c = 0; c < nchan; c++) {
			int reflayer = s->rgb_equal && nchan == 3 ? reglayer : c;
			double loc = est_loc(stats[i][c], s->lite_norm), sca = est_scale(stats[i][c], s->lite_norm);
			double pscale = 1.0, poffset = 0.0, pmul = 1.0;
			switch (s->norm) {
			case ADDITIVE_SCALING:
				pscale = sca == 0.0 ? 1.0 : s->base_scale[reflayer] / sca;
				/* no break */
			case ADDITIVE:
				poffset = pscale * loc - s->base_loc[reflayer];
				break;
			case MULTIPLICATIVE_SCALING:
				pscale = sca == 0.0 ? 1.0 : s->base_scale[reflayer] / sca;
				/* no break */
			case MULTIPLICATIVE:
				pmul = loc == 0.0 ? 1.0 : s->base_loc[reflayer] / loc;
				break;
			default:
				break;
			}
			b.pscale[c][k] = pscale;
			b.poffset[c][k] = poffset;
			b.pmul[c][k] = pmul;
			b.sframe[c][k] = stats[i][c]->bgnoise * fabs(pscale * pmul);
			noise_raw[c * nb + k] = lsf_noise_raw_weight(pscale, stats[i][c]->bgnoise);
		}
		wfwhm[k] = has_reg ? seq->regparam[reglayer][i].weighted_fwhm : 0.f;
		nstars[k] = has_reg ? seq->regparam[reglayer][i].number_of_stars : 0;
	}

	/* 4. weights, with constants frozen from the first batch */
	if ((s->weighting == WFWHM_WEIGHT || s->weighting == NBSTARS_WEIGHT) && !has_reg) {
		siril_log_error(_("Sequence does not have registration info, cannot use weighing by %s, aborting\n"),
				lsf_weight_name(s->weighting));
		free(noise_raw); free(wfwhm); free(nstars); free(cand);
		goto end;
	}
	if (args->init && lsf_freeze_weights(s, nb, nchan, noise_raw, wfwhm, nstars)) {
		free(noise_raw); free(wfwhm); free(nstars); free(cand);
		goto end;
	}
	int kept = 0;
	for (int k = 0; k < nb; k++) {
		gboolean zero = FALSE;
		double w[3];
		for (int c = 0; c < nchan; c++) {
			w[c] = lsf_frame_weight(s, c, noise_raw[c * nb + k], wfwhm[k], nstars[k]);
			zero |= !(w[c] > 0.0);
		}
		if (zero) {
			siril_log_warning(_("LSF: frame %d has a weight of 0 and is not merged\n"),
					seq->imgparam[cand[k]].filenum);
			continue;
		}
		b.idx[kept] = cand[k];
		for (int c = 0; c < nchan; c++) {
			b.weight[c][kept] = w[c];
			b.pscale[c][kept] = b.pscale[c][k];
			b.poffset[c][kept] = b.poffset[c][k];
			b.pmul[c][kept] = b.pmul[c][k];
			b.sframe[c][kept] = b.sframe[c][k];
		}
		kept++;
	}
	free(noise_raw); free(wfwhm); free(nstars); free(cand);
	b.nb = nb = kept;
	if (nb == 0) {
		siril_log_error(_("LSF: no frame to merge in this batch\n"));
		goto end;
	}
	if (s->rej == GESDT && nb < 3)
		siril_log_warning(_("LSF: the Generalized Extreme Studentized Deviate Test needs at least three frames, "
					"only the history test applies to this batch\n"));
	if (args->init && nb < LSF_LOW_STACK_WARN)
		siril_log_warning(_("LSF: only %d frames in the first batch of this new LSF file (fewer than %d). "
					"Per-pixel scatter estimated from so few frames is unreliable, and any rejection "
					"decision made on it stays in the stack.\n"), nb, LSF_LOW_STACK_WARN);

	/* 5. open the images (gives dates and the reference metadata too) */
	sargs.seq = seq;
	sargs.nb_images_to_stack = nb;
	sargs.image_indices = b.idx;
	sargs.ref_image = b.idx[0];
	for (int k = 0; args->init && k < nb; k++)
		if (b.idx[k] == args->ref_image)
			sargs.ref_image = args->ref_image;
	sargs.reglayer = -1;
	sargs.weighting_type = NO_WEIGHT;
	int bitpix = 0, naxis = 0;
	long naxes[3] = { 0, 0, 1 };
	if (stack_open_all_files(&sargs, &bitpix, &naxis, naxes, &list_date, &refmeta)) {
		siril_log_error(_("LSF: could not open the frames of the batch\n"));
		goto end;
	}
	opened = TRUE;
	if (sargs.drizzle) {
		siril_log_error(_("LSF: drizzle sequences are not supported yet\n"));
		goto end;
	}
	naxes[2] = nchan;
	gboolean ushort_input = get_data_type(bitpix) == DATA_USHORT;

	rargs.type_of_rejection = (s->rej == GESDT && nb < 3) ? NO_REJEC : s->rej;
	rargs.sig[0] = (float) s->sig[0];
	rargs.sig[1] = (float) s->sig[1];
	if (rargs.type_of_rejection == GESDT) {
		int max_outliers = (int) floor(nb * rargs.sig[0]);
		rargs.critical_value = malloc(MAX(max_outliers, 1) * sizeof(float));
		for (int j = 0, size = nb; j < max_outliers; j++, size--) {
			float t_dist = gsl_cdf_tdist_Pinv(1 - rargs.sig[1] / (2 * size), size - 2);
			float numerator = (size - 1) * t_dist;
			float denominator = sqrtf(size) * sqrtf(size - 2 + (t_dist * t_dist));
			rargs.critical_value[j] = numerator / denominator;
		}
	}

	/* history test constants (Section 10.7) */
	double kl = s->hist_low, kh = s->hist_high;
	if (lsf_rejection_has_sigma(s->rej)) {
		kl = MIN(kl, s->sig[0]);
		kh = MIN(kh, s->sig[1]);
	}
	double clip_c = lsf_clip_correction(kl, kh);
	double margin = MAX(s->hist_margin, 1.0);

	/* 6. blocks of rows, memory-bounded */
#ifdef _OPENMP
	nb_threads = com.max_thread;
	if (nb_threads > 1 && (seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ) && !fits_is_reentrant())
		nb_threads = 1;
#endif
	guint64 mem = (guint64) get_max_memory_in_MB() * BYTES_IN_A_MB;
	guint64 row_bytes = (guint64) nb * naxes[0] * (sizeof(float) + (ushort_input ? sizeof(WORD) : 0));
	long max_rows = (long) (mem / (MAX(row_bytes, 1) * nb_threads));
	if (max_rows < 1) {
		siril_log_error(_("LSF: not enough memory to process the batch\n"));
		goto end;
	}
	if (max_rows > naxes[1]) max_rows = naxes[1];
	long largest_block_height = 0;
	int nb_blocks = 0;
	if (stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block_height, &nb_blocks))
		goto end;
	size_t block_pix = (size_t) largest_block_height * naxes[0];
	tdata = calloc(nb_threads, sizeof(lsf_thread_data));
	for (int t = 0; t < nb_threads; t++) {
		if (alloc_thread_data(&tdata[t], nb, block_pix, ushort_input, rargs.type_of_rejection)) {
			PRINT_ALLOC_ERR;
			goto end;
		}
	}

	gui_iface.set_progress(PROGRESS_RESET, _("LSF: merging the batch"));
	guint64 irej[3][2] = { { 0 } }, hrej[3][2] = { { 0 } };
	int err = 0, done_blocks = 0;
	gboolean saturated = FALSE;
	const long rx = naxes[0], ry = naxes[1];

#ifdef _OPENMP
#pragma omp parallel for num_threads(nb_threads) schedule(dynamic) if (nb_threads > 1)
#endif
	for (int bi = 0; bi < nb_blocks; bi++) {
		if (err || !processing_should_continue()) {
			err = 1;
			continue;
		}
		int tid = 0;
#ifdef _OPENMP
		tid = omp_get_thread_num();
#endif
		lsf_thread_data *td = &tdata[tid];
		struct _image_block *blk = blocks + bi;
		int layer = (int) blk->channel;
		td->db.layer = layer;
		guint64 brej[2] = { 0, 0 }, bhrej[2] = { 0, 0 };
		gboolean bsat = FALSE;

		/* read the block from all frames */
		rectangle area = { 0, (int) blk->start_row, (int) rx, (int) blk->height };
		size_t npix = (size_t) blk->height * rx;
		for (int f = 0; f < nb && !err; f++) {
			void *buffer = ushort_input ? (void *) td->wbuf : (void *) td->pix[f];
			if (seq_opened_read_region(seq, layer, b.idx[f], buffer, &area, tid)) {
				siril_log_error(_("Error reading one of the image areas (%d: %d %d %d %d)\n"),
						b.idx[f] + 1, area.x, area.y, area.w, area.h);
				err = 1;
				break;
			}
			if (ushort_input)
				for (size_t p = 0; p < npix; p++)
					td->pix[f][p] = (float) (td->wbuf[p] * INV_USHRT_MAX_DOUBLE);
		}
		if (err) continue;

		float *stack = (float *) td->db.stack;
		float *o_stack = (float *) td->db.o_stack;
		float *mean_h = h->plane[LSF_PL_MEAN][layer], *w_h = h->plane[LSF_PL_W][layer];
		float *m2_h = h->plane[LSF_PL_M2][layer], *w2_h = h->plane[LSF_PL_W2][layer];
		WORD *count_h = h->count[layer];

		for (long y = 0; y < blk->height; y++) {
			/* line y of the block is stored at ry - y - 1, as in
			 * stack_mean_with_rejection(), so the Stack matches Siril's */
			size_t pdata_idx = (size_t) (ry - (blk->start_row + y) - 1) * rx;
			size_t line_idx = (size_t) y * rx;
			for (long x = 0; x < rx; x++, pdata_idx++) {
				for (int f = 0; f < nb; f++) {
					float v = td->pix[f][line_idx + x];
					stack[f] = (float) normalize_sample(s->norm, v, b.pscale[layer][f],
							b.poffset[layer][f], b.pmul[layer][f]);
				}
				int rej[2] = { 0, 0 };
				int kept_pixels = apply_rejection_float(&td->db, nb, &rargs, rej);
				brej[0] += rej[0];
				brej[1] += rej[1];
				if (kept_pixels <= 0)
					continue;	// no data here: nothing is added
				/* survivors, identified as mean_and_reject() does: its weighted
				 * branch takes the values of o_stack within the range of the kept
				 * samples, its unweighted branch the kept samples themselves. They
				 * differ only when the stop rule of the rejection (N - r <= 4) left
				 * a gap in the kept range, then the kept samples are matched back
				 * to their frames. Sigmedian keeps the range test: its kept samples
				 * include the medians that replaced the outliers. */
				float pmin = FLT_MAX, pmax = -FLT_MAX;
				for (int f = 0; f < kept_pixels; f++) {
					if (pmin > stack[f]) pmin = stack[f];
					if (pmax < stack[f]) pmax = stack[f];
				}
				gboolean *sel = td->sel;
				int nsel = 0;
				for (int f = 0; f < nb; f++) {
					float val = o_stack[f];
					sel[f] = val >= pmin && val <= pmax && val != 0.f;
					nsel += sel[f];
				}
				if (nsel != kept_pixels && s->weighting == NO_WEIGHT && s->rej != SIGMEDIAN) {
					memset(sel, 0, nb * sizeof(gboolean));
					for (int k = 0; k < kept_pixels; k++)
						for (int f = 0; f < nb; f++)
							if (!sel[f] && o_stack[f] == stack[k]) {
								sel[f] = TRUE;
								break;
							}
				}
				double Wh = w_h[pdata_idx];
				int n = 0;
				for (int f = 0; f < nb; f++) {
					float val = o_stack[f];
					if (!sel[f])
						continue;
					double wf = b.weight[layer][f];
					if (Wh > 0.0) {	// history test, Section 10.7
						double sp2 = lsf_sigma_pred2(s->weighting, Wh, m2_h[pdata_idx],
								w2_h[pdata_idx], count_h[pdata_idx], wf);
						double sh = lsf_sigma_hist(sp2, clip_c, b.sframe[layer][f], margin);
						if (sh > 0.0) {
							double z = ((double) val - mean_h[pdata_idx]) / sh;
							if (z < -s->hist_low) { bhrej[0]++; continue; }
							if (z > s->hist_high) { bhrej[1]++; continue; }
						}
					}
					td->xbuf[n] = val;
					td->wbufd[n] = wf;
					n++;
				}
				if (n == 0)
					continue;
				double mb, wb, m2b, w2b;
				lsf_batch_stats(td->xbuf, td->wbufd, n, &mb, &wb, &m2b, &w2b);
				lsf_merge(&mean_h[pdata_idx], &w_h[pdata_idx], &m2_h[pdata_idx], &w2_h[pdata_idx],
						&count_h[pdata_idx], mb, wb, m2b, w2b, n, &bsat);
			}
		}
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			irej[layer][0] += brej[0];
			irej[layer][1] += brej[1];
			hrej[layer][0] += bhrej[0];
			hrej[layer][1] += bhrej[1];
			saturated |= bsat;
			done_blocks++;
			gui_iface.set_progress((double) done_blocks / nb_blocks, NULL);
		}
	}
	if (err) {
		siril_log_error(_("LSF: merging the batch failed or was cancelled\n"));
		goto end;
	}
	h->count_saturated |= saturated;

	double nb_tot = (double) rx * ry * nb;
	for (int c = 0; c < nchan; c++) {
		siril_log_message(_("LSF pixel rejection in channel #%d: batch %.3lf%% - %.3lf%%, history %.3lf%% - %.3lf%%\n"),
				c, irej[c][0] / nb_tot * 100.0, irej[c][1] / nb_tot * 100.0,
				hrej[c][0] / nb_tot * 100.0, hrej[c][1] / nb_tot * 100.0);
		args->rejected[c][0] = irej[c][0];
		args->rejected[c][1] = irej[c][1];
		args->hist_rejected[c][0] = hrej[c][0];
		args->hist_rejected[c][1] = hrej[c][1];
	}

	/* 7. ledger rows and the header of the Stack */
	for (int k = 0; k < nb; k++) {
		lsf_ledger_row *r = g_new0(lsf_ledger_row, 1);
		fill_ledger_row(r, args, &b, k, &metas[b.idx[k]], stats[b.idx[k]]);
		g_ptr_array_add(args->rows, r);
	}
	if (args->init) {
		copy_fits_metadata(&refmeta, args->primary);
		args->primary->keywords.livetime = 0.0;
		args->primary->keywords.expstart = args->primary->keywords.expend = DEFAULT_DOUBLE_VALUE;
		if (args->primary->keywords.date_obs) {
			g_date_time_unref(args->primary->keywords.date_obs);
			args->primary->keywords.date_obs = NULL;
		}
		g_slist_free_full(args->primary->history, g_free);
		args->primary->history = NULL;
	}
	fits *p = args->primary;
	double livetime = dbl_known(p->keywords.livetime) ? p->keywords.livetime : 0.0;
	for (int k = 0; k < nb; k++) {
		const fits *m = &metas[b.idx[k]];
		if (dbl_known(m->keywords.exposure))
			livetime += m->keywords.exposure;
		reconcile_first_frame_keywords(p, m);
	}
	p->keywords.livetime = livetime;
	p->keywords.stackcnt = (guint) (args->prev_frames + nb);
	if (list_date) {
		fits tmp = { 0 };
		compute_date_time_keywords(list_date, &tmp);	// sorts the list
		list_date = g_list_first(list_date);
		if (!dbl_known(p->keywords.expstart) || tmp.keywords.expstart < p->keywords.expstart)
			p->keywords.expstart = tmp.keywords.expstart;
		if (!dbl_known(p->keywords.expend) || tmp.keywords.expend > p->keywords.expend)
			p->keywords.expend = tmp.keywords.expend;
		GDateTime *first = ((DateEvent *) list_date->data)->date_obs;
		if (!p->keywords.date_obs || g_date_time_compare(first, p->keywords.date_obs) < 0) {
			if (p->keywords.date_obs)
				g_date_time_unref(p->keywords.date_obs);
			p->keywords.date_obs = g_date_time_ref(first);
		}
	}
	/* the Stack is already RGB data and has no display cutoffs (Sections 4.1, 7) */
	p->keywords.bayer_pattern[0] = '\0';
	p->keywords.lo = p->keywords.hi = 0;
	p->keywords.flo = p->keywords.fhi = 0.f;
	if (!p->keywords.row_order[0])
		g_strlcpy(p->keywords.row_order, "BOTTOM-UP", sizeof p->keywords.row_order);
	g_strlcpy(p->keywords.program, "Siril " PACKAGE_VERSION, sizeof p->keywords.program);
	p->history = g_slist_append(p->history, g_strdup_printf("LSF generation %d: batch %d, %d frames from %s",
				args->gen, args->batch, nb, args->seq_name));
	/* rejected samples in percent, low and high: batch-local then history */
	for (int c = 0; c < nchan; c++)
		p->history = g_slist_append(p->history, g_strdup_printf("LSF gen %d ch %d rejected %%: %.4f %.4f hist %.4f %.4f",
					args->gen, c, irej[c][0] / nb_tot * 100.0, irej[c][1] / nb_tot * 100.0,
					hrej[c][0] / nb_tot * 100.0, hrej[c][1] / nb_tot * 100.0));
	retval = 0;

end:
	if (opened)
		for (int k = 0; k < b.nb; k++)
			seq_close_image(seq, b.idx[k]);
	if (tdata) {
		for (int t = 0; t < nb_threads; t++)
			free_thread_data(&tdata[t]);
		free(tdata);
	}
	free(blocks);
	free(rargs.critical_value);
	g_list_free_full(list_date, (GDestroyNotify) free_list_date);
	clearfits(&refmeta);
	if (metas) {
		for (int i = 0; i < seq->number; i++)
			clearfits(&metas[i]);
		free(metas);
	}
	if (stats) {
		for (int i = 0; i < seq->number; i++) {
			if (!stats[i]) continue;
			for (int c = 0; c < 3; c++)
				if (stats[i][c]) free_stats(stats[i][c]);
			free(stats[i]);
		}
		free(stats);
	}
	free_batch_frames(&b);
	gui_iface.set_progress(retval ? PROGRESS_RESET : PROGRESS_DONE, NULL);
	return retval;
}

/* HDU 0 from MEAN: output-normalized when enabled, otherwise clamped to [0,1]
 * like the 32-bit output of stack_mean_with_rejection(), 0 where W = 0
 * (Sections 3, 4.1). MEAN itself is not clamped. */
int lsf_build_primary(fits *primary, lsf_settings *s, const lsf_history *h) {
	size_t n = (size_t) h->rx * h->ry;
	float *data = malloc(n * h->nchan * sizeof(float));
	if (!data) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	double mini = 0.0, maxi = 1.0;
	gboolean norm = s->outnorm;
	if (norm && !lsf_outnorm_range(h->plane[LSF_PL_MEAN], h->plane[LSF_PL_W], n, h->nchan, &mini, &maxi)) {
		siril_log_warning(_("LSF: output normalization skipped, no valid dynamic range\n"));
		norm = FALSE;
	}
	double range = maxi - mini;
	for (int c = 0; c < h->nchan; c++) {
		const float *m = h->plane[LSF_PL_MEAN][c], *w = h->plane[LSF_PL_W][c];
		float *out = data + c * n;
		for (size_t i = 0; i < n; i++)
			out[i] = lsf_stack_value(m[i], w[i], norm, mini, range);
	}
	s->outnorm_min = norm ? mini : 0.0;
	s->outnorm_max = norm ? maxi : 0.0;
	s->outnorm = norm;
	free(primary->data);
	free(primary->fdata);
	primary->data = NULL;
	for (int c = 0; c < 3; c++)
		primary->pdata[c] = NULL;
	primary->fdata = data;
	for (int c = 0; c < 3; c++)
		primary->fpdata[c] = data + (c < h->nchan ? c : 0) * n;
	primary->type = DATA_FLOAT;
	primary->bitpix = FLOAT_IMG;
	primary->orig_bitpix = FLOAT_IMG;
	primary->rx = h->rx;
	primary->ry = h->ry;
	primary->naxes[0] = h->rx;
	primary->naxes[1] = h->ry;
	primary->naxes[2] = h->nchan;
	primary->naxis = h->nchan == 3 ? 3 : 2;
	return 0;
}
