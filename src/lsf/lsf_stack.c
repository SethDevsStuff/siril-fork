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

/* LSF batch engine (FITS profile v0.4, Sections 10.2, 10.3, 10.6).
 *
 * The batch, a sequence already on the LSF grid, is read in blocks of rows
 * like stack_mean_with_rejection() does. Each pixel goes through Siril's own
 * normalization formulas (with coefficients computed against the frozen
 * baseline of the file), Siril's batch-local rejection
 * (apply_rejection_float), then the LSF history test. The accepted samples
 * are merged into the accumulator with the weighted pairwise update.
 * Siril's stacking code is not modified. */

#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_cdf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "algos/star_finder.h"
#include "algos/statistics.h"
#include "core/gui_iface.h"
#include "core/OS_utils.h"
#include "core/processing.h"
#include "core/proto.h"
#include "core/siril.h"
#include "core/siril_alloc.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "io/sequence.h"
#include "lsf/lsf.h"
#include "lsf/lsf_fits.h"
#include "lsf/lsf_math.h"
#include "lsf/lsf_settings.h"
#include "lsf/lsf_stack.h"
#include "lsf/lsf_wcs.h"
#include "lsf/lsf_weights.h"
#include "opencv/opencv.h"
#include "registration/registration.h"

void lsf_ledger_row_free(gpointer p) {
	lsf_ledger_row *r = p;
	if (!r) return;
	g_free(r->frame_id);
	g_free(r->filename);
	g_free(r->date_obs);
	g_free(r->filter);
	g_free(r->orighdr);
	g_free(r->align_json);
	g_free(r->reg_dist);
	g_free(r);
}

/* per-frame data of the batch */
typedef struct {
	int nb;			// frames kept
	int *idx;		// sequence index of each kept frame
	double *pscale[3], *poffset[3], *pmul[3];
	double *weight[3];
	double *sframe[3];	// frame noise in accumulator units (sigma_frame)
	/* star measures on the grid against LSF_WREF (Section 10.3) */
	float *fwhm, *wfwhm, *round;
	int *nstars;
	/* the samples that count as data, lsf_edge_limits() of each frame:
	 * [frame * ry + y] and [frame * rx + x] (Section 10.5) */
	int *row_lo, *row_hi, *col_lo, *col_hi;
} batch_frames;

static void free_batch_frames(batch_frames *b) {
	free(b->row_lo); free(b->row_hi); free(b->col_lo); free(b->col_hi);
	free(b->fwhm);
	free(b->wfwhm);
	free(b->round);
	free(b->nstars);
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

/* layer holding the registration data of the batch, as the stack command
 * finds it; green or mono when there is none */
static int batch_reg_layer(const sequence *seq) {
	int layer = get_registration_layer(seq);
	return layer >= 0 ? layer : (seq->nb_layers == 3 ? 1 : 0);
}

static gboolean differ_dbl(double a, double b) {
	return dbl_known(a) && dbl_known(b) && fabs(a - b) > 1e-6 * MAX(fabs(a), fabs(b));
}

static gboolean differ_str(const char *a, const char *b) {
	return a[0] != '\0' && b[0] != '\0' && strcmp(a, b);
}

/* Section 7, keywords set at initialization then checked on append: the
 * first that differs, or NULL. The GUI checks a batch with it beforehand. */
const char *lsf_frame_mismatch(const fits *primary, const fits *frame) {
	const fkeywords *p = &primary->keywords, *f = &frame->keywords;
	if (p->binning_x && f->binning_x && (p->binning_x != f->binning_x || p->binning_y != f->binning_y))
		return "XBINNING/YBINNING";
	if (differ_dbl(p->pixel_size_x, f->pixel_size_x) || differ_dbl(p->pixel_size_y, f->pixel_size_y))
		return "XPIXSZ/YPIXSZ";
	/* plate solving replaces FOCALLEN with the focal length it measures, a
	 * few percent from the nominal one; a reducer or another telescope
	 * differs by much more */
	if (dbl_known(p->focal_length) && dbl_known(f->focal_length)
			&& fabs(p->focal_length - f->focal_length) > LSF_FOCAL_TOLERANCE * MAX(p->focal_length, f->focal_length))
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
	int reglayer = batch_reg_layer(seq);
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
	/* the distortion model, as ALIGN names it and as REG_DIST records the
	 * SIP terms applied (Section 9) */
	gchar *disto = g_strdup("");
	if (args->src_seq && args->disto && args->disto->index != DISTO_UNDEF) {
		g_free(disto);
		disto = g_strdup_printf(",\"disto\":\"%s\"", args->disto->index == DISTO_MASTER ? "master" : "file");
		r->reg_dist = lsf_disto_cards(args->disto, (fits *) meta);
	}
	if (!args->src_seq)
		r->align_json = g_strdup("{\"reference\":\"none\"}");
	else if (args->xfrm_base)
		r->align_json = g_strdup_printf("{\"reference\":\"sequence\",\"framing\":\"%s\",\"transform\":\"%s\",\"interpolation\":%d%s}",
				args->framing ? args->framing : "current", transf[t], args->settings->interpolation, disto);
	else r->align_json = g_strdup_printf("{\"reference\":\"stack\",\"method\":\"%s\",\"framing\":\"%s\",\"transform\":\"%s\","
				"\"interpolation\":%d,\"min_pairs\":%d%s}", args->reg_method ? args->reg_method : "stars",
				args->framing ? args->framing : "current", transf[t < 0 ? 3 : t], args->settings->interpolation,
				args->settings->min_pairs, disto);
	g_free(disto);

	regdata *q = reg ? reg : src_reg;
	r->regq = q ? (float) q->quality : NAN;
	r->regbg = q ? q->background_lvl : NAN;
	/* the measures the weights used, on the grid (Section 10.3) */
	r->fwhm = b->fwhm[k];
	r->wfwhm = b->wfwhm[k];
	r->round = b->round[k];
	r->nstars = b->nstars[k];

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
	t->mem = siril_malloc((size_t) nb * block_pix * sizeof(float));
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
	free(t->pix); siril_free(t->mem); free(t->wbuf); free(t->xbuf); free(t->wbufd); free(t->sel);
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

/***** star measures against LSF_WREF (Section 10.3) *****/

/* the radius of a star pair: the RANSAC threshold of Siril's star alignment
 * (defaultRANSACReprojThreshold in opencv.cpp), within which it keeps a pair */
#define LSF_PAIR_RADIUS 3.0

/* the reference stars where the frame, on the grid, has data: the others
 * cannot pair, and a frame framed differently would be counted as worse */
static int count_covered(const fits *fit, int layer, const lsf_wref *ref) {
	int covered = 0;
	size_t plane = (size_t) fit->rx * fit->ry * (fit->naxes[2] > 1 ? layer : 0);
	for (int j = 0; j < ref->n; j++) {
		/* FITS pixels to the pixel holding them, in FITS order */
		long px = lround(ref->x[j] - 1.0), py = lround(ref->y[j] - 1.0);
		if (px < 0 || py < 0 || px >= fit->rx || py >= fit->ry)
			continue;
		size_t at = plane + (size_t) py * fit->rx + px;
		double v = fit->type == DATA_FLOAT ? fit->fdata[at] : fit->data[at];
		if (v != 0.0 && isfinite(v))
			covered++;
	}
	return covered;
}

/* lsf_edge_limits() of a frame on the grid, a pixel being data when it is
 * not 0 in some channel. Returns 0, or 1 on error. */
static int frame_edge_limits(const fits *fit, int margin, int *row_lo, int *row_hi, int *col_lo, int *col_hi) {
	size_t n = (size_t) fit->rx * fit->ry;
	guint8 *has = calloc(n, 1);
	if (!has) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	for (int c = 0; c < (int) fit->naxes[2]; c++) {
		for (size_t p = 0; p < n; p++)
			has[p] |= fit->type == DATA_FLOAT ? fit->fdata[c * n + p] != 0.f : fit->data[c * n + p] != 0;
	}
	lsf_edge_limits(has, fit->rx, fit->ry, margin, row_lo, row_hi, col_lo, col_hi);
	free(has);
	return 0;
}

/* the stars of frame i of seq, detected with sf on layer, in Siril
 * coordinates (FITS order, 0-based), and its FWHM, roundness, and the number
 * of reference stars where it has data. With edges, the edge limits of the
 * frame for its batch slot k too. Returns the number of stars, 0 when none,
 * -1 on error. */
static int detect_frame_stars(sequence *seq, int i, struct starfinder_data *sf, int threads, const lsf_wref *ref,
		batch_frames *edges, int k, int margin, double **x, double **y, float *fwhm, float *round, int *covered) {
	fits fit = { 0 };
	*x = *y = NULL;
	if (seq_read_frame(seq, i, &fit, FALSE, -1)) {
		siril_log_error(_("LSF: could not read frame %d\n"), seq->imgparam[i].filenum);
		return -1;
	}
	if (edges && frame_edge_limits(&fit, margin, edges->row_lo + (size_t) k * fit.ry, edges->row_hi + (size_t) k * fit.ry,
				edges->col_lo + (size_t) k * fit.rx, edges->col_hi + (size_t) k * fit.rx)) {
		clearfits(&fit);
		return -1;
	}
	*covered = ref ? count_covered(&fit, sf->layer, ref) : 0;
	psf_star **stars = NULL;
	int nb = 0;
	struct starfinder_data *found = findstar_image_worker(sf, -1, i, &fit, NULL, threads);
	if (found) {
		stars = *found->stars;
		nb = *found->nb_stars;
		free(found);
	}
	*fwhm = *round = NAN;
	if (stars && nb > 0) {
		float FWHMx, FWHMy, B;
		char *units;
		FWHM_stats(stars, nb, seq->bitpix, &FWHMx, &FWHMy, &units, &B, NULL, 0.);
		*fwhm = FWHMx;
		*round = FWHMy / FWHMx;
		*x = malloc(nb * sizeof(double));
		*y = malloc(nb * sizeof(double));
		if (!*x || !*y) {
			PRINT_ALLOC_ERR;
			nb = -1;
		}
		for (int k = 0; k < nb; k++)
			display_to_siril(stars[k]->xpos, stars[k]->ypos, *x + k, *y + k, fit.ry);
	} else nb = 0;
	if (stars)
		free_fitted_stars(stars);
	clearfits(&fit);
	if (nb < 0) {
		free(*x);
		free(*y);
		*x = *y = NULL;
	}
	return nb;
}

/* pairs between a frame's stars and the reference stars, both on the grid,
 * each reference star used once, within LSF_PAIR_RADIUS */
static int count_pairs(const double *x, const double *y, int nb, const lsf_wref *ref) {
	gboolean *used = calloc(ref->n, sizeof(gboolean));
	if (!used) {
		PRINT_ALLOC_ERR;
		return 0;
	}
	int pairs = 0;
	const double r2 = LSF_PAIR_RADIUS * LSF_PAIR_RADIUS;
	for (int i = 0; i < nb; i++) {
		int best = -1;
		double bd = r2;
		for (int j = 0; j < ref->n; j++) {
			if (used[j])
				continue;
			/* LSF_WREF holds FITS pixels: Siril's coordinates + 0.5 */
			double dx = ref->x[j] - 0.5 - x[i], dy = ref->y[j] - 0.5 - y[i], dd = dx * dx + dy * dy;
			if (dd < bd) {
				bd = dd;
				best = j;
			}
		}
		if (best >= 0) {
			used[best] = TRUE;
			pairs++;
		}
	}
	free(used);
	return pairs;
}

/* The star measures of the frames to merge, on the grid, as the weights use
 * them (Section 10.3). On initialization the stars of the reference frame
 * become LSF_WREF, frozen. Each frame's wFWHM is Siril's star alignment
 * formula with the reference stars as reference: 2 FWHM (fitted - pairs) /
 * fitted + FWHM, fitted being the reference stars where the frame has data,
 * capped as star alignment caps them. Measured on every frame the same way,
 * whatever registered it, however deep the Stack is and however it is
 * framed, so that weights stay comparable across appends. Returns 0, or 1
 * on error. The edge limits of the frames are found on the way. */
static int measure_frames(struct lsf_batch_args *args, const int *cand, int nb, batch_frames *b, int threads) {
	sequence *seq = args->seq;
	lsf_wref *ref = args->wref;
	star_finder_params sf = com.pref.starfinder_conf;
	lsf_findstar_from_json(args->settings->findstar_json, &sf);
	struct starfinder_data sfd = { 0 };
	sfd.im.from_seq = seq;
	sfd.params = &sf;
	sfd.layer = batch_reg_layer(seq);
	sfd.keep_stars = TRUE;
	sfd.max_stars_fitted = MAX_STARS_FITTED;
	b->fwhm = malloc(nb * sizeof(float));
	b->wfwhm = malloc(nb * sizeof(float));
	b->round = malloc(nb * sizeof(float));
	b->nstars = malloc(nb * sizeof(int));
	b->row_lo = malloc((size_t) nb * seq->ry * sizeof(int));
	b->row_hi = malloc((size_t) nb * seq->ry * sizeof(int));
	b->col_lo = malloc((size_t) nb * seq->rx * sizeof(int));
	b->col_hi = malloc((size_t) nb * seq->rx * sizeof(int));
	if (!b->fwhm || !b->wfwhm || !b->round || !b->nstars || !b->row_lo || !b->row_hi || !b->col_lo || !b->col_hi) {
		PRINT_ALLOC_ERR;
		return 1;
	}
	/* the frames are read here anyway: their edges are found on the way */
	const int margin = lsf_edge_margin(args->settings->interpolation);
	if (args->init && ref->n == 0) {
		double *x, *y;
		float fwhm, round;
		int covered;
		int n = detect_frame_stars(seq, args->ref_image, &sfd, threads, NULL, NULL, 0, 0, &x, &y, &fwhm, &round,
				&covered);
		if (n <= 0) {
			siril_log_error(_("LSF: no star found in the reference frame %d, the star measures of the frames cannot "
						"be referred to it\n"), seq->imgparam[args->ref_image].filenum);
			free(x);
			free(y);
			return 1;
		}
		for (int k = 0; k < n; k++) {
			x[k] += 0.5;	// FITS pixels
			y[k] += 0.5;
		}
		ref->n = n;
		ref->x = x;
		ref->y = y;
		siril_log_message(_("LSF: %d stars of the reference frame %d are the reference of the wFWHM of every frame\n"),
				n, seq->imgparam[args->ref_image].filenum);
	}
	gui_iface.set_progress(PROGRESS_RESET, _("LSF: measuring the stars of the frames"));
	for (int k = 0; k < nb; k++) {
		double *x, *y;
		int covered;
		int n = detect_frame_stars(seq, cand[k], &sfd, threads, ref, b, k, margin, &x, &y, &b->fwhm[k], &b->round[k],
				&covered);
		if (n < 0)
			return 1;
		int pairs = n > 0 ? count_pairs(x, y, n, ref) : 0;
		free(x);
		free(y);
		int fitted = MIN(covered, MAX_STARS_FITTED);
		b->nstars[k] = n;
		b->wfwhm[k] = n > 0 && fitted > 0 ? (float) (2.0 * b->fwhm[k] * MAX(fitted - pairs, 0) / fitted + b->fwhm[k]) : 0.f;
		siril_log_debug("LSF: frame %d, %d stars, %d pairs with the %d reference stars it covers, FWHM %.3f, wFWHM %.3f\n",
				seq->imgparam[cand[k]].filenum, n, pairs, fitted, b->fwhm[k], b->wfwhm[k]);
		gui_iface.set_progress((double) (k + 1) / nb, NULL);
		if (!processing_should_continue())
			return 1;
	}
	return 0;
}

int lsf_process_batch(struct lsf_batch_args *args) {
	sequence *seq = args->seq;
	lsf_settings *s = args->settings;
	lsf_history *h = args->hist;
	int nchan = seq->nb_layers;
	/* norm.rgb_equal equalizes to green: the baseline is frozen per channel,
	 * so the reference channel cannot follow the registration layer */
	int rgb_ref = 1;
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
	const int max_threads = MAX(args->max_threads, 1);	// the budget of the job, not com.max_thread
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
	/* rejection maps asked for now: they count from this generation on (Section 5.4) */
	if (s->rejmaps && !h->rej[0][0]) {
		if (lsf_history_alloc_rej(h, args->gen))
			return 1;
		if (!args->init)
			siril_log_message(_("LSF: rejection maps enabled, counting from generation %d\n"), args->gen);
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
			const char *bad = lsf_frame_mismatch(args->primary, &metas[i]);
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
		if (!stats[i]) {
			PRINT_ALLOC_ERR;
			free(cand);
			goto end;
		}
		if (compute_all_channels_statistics_seqimage(seq, i, NULL, s->lite_norm ? STATS_LITENORM : STATS_NORM,
					max_threads, -1, stats[i])) {
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
	/* caches the statistics, unless the selection was changed in memory */
	if (!args->keep_seq_file)
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
	gboolean alloc_ok = b.idx && noise_raw && wfwhm && nstars;
	for (int c = 0; c < nchan; c++)
		alloc_ok = alloc_ok && b.pscale[c] && b.poffset[c] && b.pmul[c] && b.weight[c] && b.sframe[c];
	if (!alloc_ok) {
		PRINT_ALLOC_ERR;
		free(noise_raw); free(wfwhm); free(nstars); free(cand);
		goto end;
	}
	for (int k = 0; k < nb; k++) {
		int i = cand[k];
		for (int c = 0; c < nchan; c++) {
			int reflayer = s->rgb_equal && nchan == 3 ? rgb_ref : c;
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
	}
	/* the star measures of the weights, against the frozen reference stars,
	 * whatever the weighting, so that the ledger's are comparable */
	if (measure_frames(args, cand, nb, &b, max_threads)) {
		free(noise_raw); free(wfwhm); free(nstars); free(cand);
		goto end;
	}
	for (int k = 0; k < nb; k++) {
		wfwhm[k] = b.wfwhm[k];
		nstars[k] = b.nstars[k];
	}

	/* 4. weights, with constants frozen from the first batch */
	if (args->init && lsf_freeze_weights(s, nb, nchan, noise_raw, wfwhm, nstars)) {
		free(noise_raw); free(wfwhm); free(nstars); free(cand);
		goto end;
	}
	/* a range frozen from the first batch clamps the frames of a batch that
	 * reach beyond it: they lose their differences, or fall to the floor */
	int better, worse, floored;
	if (!args->init && lsf_weight_range_check(s, nb, wfwhm, nstars, &better, &worse, &floored)
			&& (better > 0 || worse > 0)) {
		siril_log_warning(_("LSF: %d of the %d frames of this batch are beyond the weighting range frozen from the first "
					"batch: %d better than its best, which all get the top weight, and %d worse than its worst, "
					"which all get its weight. %d frames of the batch are at the weight floor. These weights no "
					"longer follow the frames' quality. Switching the file to the inverse weight model (lsf_append "
					"-weight_model=inverse) stops this for the next appends\n"), better + worse, nb, better, worse,
				floored);
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
		b.fwhm[kept] = b.fwhm[k];
		b.wfwhm[kept] = b.wfwhm[k];
		b.round[kept] = b.round[k];
		b.nstars[kept] = b.nstars[k];
		if (kept != k) {
			memcpy(b.row_lo + (size_t) kept * seq->ry, b.row_lo + (size_t) k * seq->ry, seq->ry * sizeof(int));
			memcpy(b.row_hi + (size_t) kept * seq->ry, b.row_hi + (size_t) k * seq->ry, seq->ry * sizeof(int));
			memcpy(b.col_lo + (size_t) kept * seq->rx, b.col_lo + (size_t) k * seq->rx, seq->rx * sizeof(int));
			memcpy(b.col_hi + (size_t) kept * seq->rx, b.col_hi + (size_t) k * seq->rx, seq->rx * sizeof(int));
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
		if (!rargs.critical_value) {
			PRINT_ALLOC_ERR;
			goto end;
		}
		for (int j = 0, size = nb; j < max_outliers; j++, size--) {
			float t_dist = gsl_cdf_tdist_Pinv(1 - rargs.sig[1] / (2 * size), size - 2);
			float numerator = (size - 1) * t_dist;
			float denominator = sqrtf(size) * sqrtf(size - 2 + (t_dist * t_dist));
			rargs.critical_value[j] = numerator / denominator;
		}
	}

	/* history test constants (Section 10.6) */
	double kl = s->hist_low, kh = s->hist_high;
	if (lsf_rejection_has_sigma(s->rej)) {
		kl = MIN(kl, s->sig[0]);
		kh = MIN(kh, s->sig[1]);
	}
	double clip_c = lsf_clip_correction(kl, kh);
	double margin = MAX(s->hist_margin, 1.0);

	/* 6. blocks of rows, memory-bounded */
#ifdef _OPENMP
	nb_threads = max_threads;
	if (nb_threads > 1 && (seq->type == SEQ_REGULAR || seq->type == SEQ_FITSEQ) && !fits_is_reentrant())
		nb_threads = 1;
#endif
	/* the accumulator stays in memory during the merge, so it counts against
	 * Siril's memory limit. With a limit set as a ratio of the available
	 * memory, the part already read may be counted twice, which only makes
	 * the blocks smaller. */
	guint64 hist_bytes = (guint64) h->rx * h->ry * h->nchan * (LSF_NB_PL * sizeof(float)
			+ (h->rej[0][0] ? 3 : 1) * sizeof(WORD));
	guint64 mem = (guint64) get_max_memory_in_MB() * BYTES_IN_A_MB;
	mem = mem > hist_bytes ? mem - hist_bytes : 0;
	guint64 row_bytes = (guint64) nb * naxes[0] * (sizeof(float) + (ushort_input ? sizeof(WORD) : 0));
	long max_rows = (long) (mem / (MAX(row_bytes, 1) * nb_threads));
	if (max_rows < 1) {
		siril_log_error(_("LSF: not enough memory to process the batch: the accumulator uses %.0f MB of the "
					"%d MB Siril may use, and reading the frames needs at least %.0f MB more\n"),
				(double) hist_bytes / BYTES_IN_A_MB, get_max_memory_in_MB(),
				(double) row_bytes * nb_threads / BYTES_IN_A_MB);
		goto end;
	}
	if (max_rows > naxes[1]) max_rows = naxes[1];
	long largest_block_height = 0;
	int nb_blocks = 0;
	if (stack_compute_parallel_blocks(&blocks, max_rows, naxes, nb_threads, &largest_block_height, &nb_blocks))
		goto end;
	size_t block_pix = (size_t) largest_block_height * naxes[0];
	tdata = calloc(nb_threads, sizeof(lsf_thread_data));
	if (!tdata) {
		PRINT_ALLOC_ERR;
		goto end;
	}
	for (int t = 0; t < nb_threads; t++) {
		if (alloc_thread_data(&tdata[t], nb, block_pix, ushort_input, rargs.type_of_rejection)) {
			PRINT_ALLOC_ERR;
			goto end;
		}
	}

	gui_iface.set_progress(PROGRESS_RESET, _("LSF: merging the batch"));
	guint64 irej[3][2] = { { 0 } }, hrej[3][2] = { { 0 } }, edge_drop[3] = { 0 }, thin[3] = { 0 };
	int err = 0, done_blocks = 0;
	gboolean saturated = FALSE, rej_saturated = FALSE;
	const long rx = naxes[0], ry = naxes[1];

#ifdef _OPENMP
#pragma omp parallel for num_threads(nb_threads) schedule(dynamic) if (nb_threads > 1)
#endif
	for (int bi = 0; bi < nb_blocks; bi++) {
		if (g_atomic_int_get(&err) || !processing_should_continue()) {
			g_atomic_int_set(&err, 1);
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
		guint64 brej[2] = { 0, 0 }, bhrej[2] = { 0, 0 }, bedge = 0, bthin = 0;
		gboolean bsat = FALSE, brsat = FALSE;

		/* read the block from all frames */
		rectangle area = { 0, (int) blk->start_row, (int) rx, (int) blk->height };
		gboolean block_err = FALSE;
		size_t npix = (size_t) blk->height * rx;
		for (int f = 0; f < nb && !block_err; f++) {
			void *buffer = ushort_input ? (void *) td->wbuf : (void *) td->pix[f];
			if (seq_opened_read_region(seq, layer, b.idx[f], buffer, &area, tid)) {
				siril_log_error(_("Error reading one of the image areas (%d: %d %d %d %d)\n"),
						b.idx[f] + 1, area.x, area.y, area.w, area.h);
				block_err = TRUE;
				break;
			}
			if (ushort_input)
				for (size_t p = 0; p < npix; p++)
					td->pix[f][p] = (float) (td->wbuf[p] * INV_USHRT_MAX_DOUBLE);
		}
		if (block_err) {
			g_atomic_int_set(&err, 1);
			continue;
		}

		float *stack = (float *) td->db.stack;
		float *o_stack = (float *) td->db.o_stack;
		float *mean_h = h->plane[LSF_PL_MEAN][layer], *w_h = h->plane[LSF_PL_W][layer];
		float *m2_h = h->plane[LSF_PL_M2][layer], *w2_h = h->plane[LSF_PL_W2][layer];
		WORD *count_h = h->count[layer];
		WORD *rlo_h = h->rej[0][layer], *rhi_h = h->rej[1][layer];	// NULL without LSF_REJ

		for (long y = 0; y < blk->height; y++) {
			/* line y of the block is stored at ry - y - 1, as in
			 * stack_mean_with_rejection(), so the Stack matches Siril's */
			const long srow = ry - (blk->start_row + y) - 1;
			size_t pdata_idx = (size_t) srow * rx;
			size_t line_idx = (size_t) y * rx;
			for (long x = 0; x < rx; x++, pdata_idx++) {
				for (int f = 0; f < nb; f++) {
					float v = td->pix[f][line_idx + x];
					/* frame-edge samples are not data (Section 10.5) */
					size_t r = (size_t) f * ry + srow, c = (size_t) f * rx + x;
					if (v != 0.f && (x < b.row_lo[r] || x > b.row_hi[r] || srow < b.col_lo[c] || srow > b.col_hi[c])) {
						v = 0.f;
						bedge++;
					}
					stack[f] = (float) normalize_sample(s->norm, v, b.pscale[layer][f],
							b.poffset[layer][f], b.pmul[layer][f]);
				}
				int rej[2] = { 0, 0 };
				int kept_pixels = apply_rejection_float(&td->db, nb, &rargs, rej);
				brej[0] += rej[0];
				brej[1] += rej[1];
				if (kept_pixels <= 0) {
					if (rlo_h) {
						lsf_count_add(&rlo_h[pdata_idx], rej[0], &brsat);
						lsf_count_add(&rhi_h[pdata_idx], rej[1], &brsat);
					}
					continue;	// no data here: nothing is added
				}
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
				gboolean testable = lsf_history_testable(Wh, w2_h[pdata_idx]);
				int n = 0, hlo = 0, hhi = 0;
				for (int f = 0; f < nb; f++) {
					float val = o_stack[f];
					if (!sel[f])
						continue;
					double wf = b.weight[layer][f];
					if (testable) {	// history test, Section 10.6
						double sp2 = lsf_sigma_pred2(s->weighting, Wh, m2_h[pdata_idx],
								w2_h[pdata_idx], count_h[pdata_idx], wf);
						double sh = lsf_sigma_hist(sp2, clip_c, b.sframe[layer][f], margin);
						if (sh > 0.0) {
							double z = ((double) val - mean_h[pdata_idx]) / sh;
							if (z < -s->hist_low) { hlo++; continue; }
							if (z > s->hist_high) { hhi++; continue; }
						}
					}
					td->xbuf[n] = val;
					td->wbufd[n] = wf;
					n++;
				}
				if (!testable && Wh > 0.0)
					bthin += n;
				bhrej[0] += hlo;
				bhrej[1] += hhi;
				/* both tests count in the rejection maps (Section 5.4) */
				if (rlo_h) {
					lsf_count_add(&rlo_h[pdata_idx], rej[0] + hlo, &brsat);
					lsf_count_add(&rhi_h[pdata_idx], rej[1] + hhi, &brsat);
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
			edge_drop[layer] += bedge;
			thin[layer] += bthin;
			saturated |= bsat;
			rej_saturated |= brsat;
			done_blocks++;
			gui_iface.set_progress((double) done_blocks / nb_blocks, NULL);
		}
	}
	if (err) {
		siril_log_error(_("LSF: merging the batch failed or was cancelled\n"));
		goto end;
	}
	h->count_saturated |= saturated;
	h->rej_saturated |= rej_saturated;

	double nb_tot = (double) rx * ry * nb;
	siril_log_message(_("LSF: %" G_GUINT64_FORMAT " frame-edge samples of channel #0 not merged (%d px at each edge "
				"of a frame inside the grid)\n"), edge_drop[0], lsf_edge_margin(s->interpolation));
	if (!args->init)
		siril_log_message(_("LSF: %" G_GUINT64_FORMAT " samples of channel #0 merged without the history test, "
					"where the history holds fewer than %g effective samples\n"), thin[0], LSF_HIST_MIN_NEFF);
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
		siril_log_warning(_("LSF: the Stack has no valid dynamic range for output normalization. "
					"Output normalization is now off for this file, and later appends will not turn it back on\n"));
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
	primary->rx = h->rx;
	primary->ry = h->ry;
	primary->naxes[0] = h->rx;
	primary->naxes[1] = h->ry;
	primary->naxes[2] = h->nchan;
	primary->naxis = h->nchan == 3 ? 3 : 2;
	/* fit_replace_buffer() frees the buffer of the other type only */
	free(primary->fdata);
	primary->fdata = NULL;
	fit_replace_buffer(primary, data, DATA_FLOAT);
	return 0;
}
