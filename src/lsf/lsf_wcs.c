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

/* LSF WCS and distortion (FITS profile v0.4, Sections 8, 9 and 10.4).
 * Siril's own WCS functions (siril_wcs.c), star detection (star_finder.c)
 * and homography fit (opencv.cpp) are used as they are. */

#include <math.h>
#include <string.h>

#include "algos/astrometry_solver.h"
#include "algos/siril_wcs.h"
#include "algos/star_finder.h"
#include "core/processing.h"
#include "core/sequence_filtering.h"
#include "core/siril.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "io/path_parse.h"
#include "io/sequence.h"
#include "lsf/lsf.h"
#include "lsf/lsf_wcs.h"
#include "opencv/opencv.h"
#include "registration/registration.h"

/* points per axis of the grid a transform is fitted on (lsf_wcs_homography) */
#define LSF_WCS_GRID 8

/* the size an image or a Siril .wcs file (NAXIS = 0, IMAGEW and IMAGEH)
 * describes, and the LSFUID of the lsf_export provenance card, if any */
static int read_import_header(const char *file, int *rx, int *ry, gchar **export_uid) {
	fitsfile *fptr = NULL;
	int status = 0, naxis = 0, nkeys = 0;
	*rx = *ry = 0;
	*export_uid = NULL;
	if (fits_open_diskfile(&fptr, file, READONLY, &status))
		return 1;
	fits_get_img_dim(fptr, &naxis, &status);
	if (!status && naxis >= 2) {
		long n[3] = { 0 };
		fits_get_img_size(fptr, 2, n, &status);
		*rx = (int) n[0];
		*ry = (int) n[1];
	} else {
		status = 0;
		fits_read_key(fptr, TINT, "IMAGEW", rx, NULL, &status);
		status = 0;
		fits_read_key(fptr, TINT, "IMAGEH", ry, NULL, &status);
	}
	status = 0;
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);
	static const char *prov = "HISTORY Exported from LSF file ";
	for (int k = 1; k <= nkeys && !status && !*export_uid; k++) {
		char card[FLEN_CARD];
		fits_read_record(fptr, k, card, &status);
		if (!status && g_str_has_prefix(card, prov)) {
			const char *u = card + strlen(prov);
			const char *end = strchr(u, ',');
			*export_uid = end ? g_strndup(u, end - u) : g_strstrip(g_strdup(u));
		}
	}
	status = 0;
	fits_close_file(fptr, &status);
	return 0;
}

/* Reads the WCS that lsf_util -wcs imports (Section 8): a plate-solved image
 * of the Stack, or a .wcs file, of the grid's size. out receives its
 * metadata, with the WCS. Returns 0, or 1 with why. */
int lsf_wcs_read_import(const char *file, const lsf_file *f, fits *out, gchar **why) {
	int rx, ry;
	gchar *uid = NULL;
	*why = NULL;
	memset(out, 0, sizeof(fits));
	if (read_import_header(file, &rx, &ry, &uid)) {
		*why = g_strdup_printf(_("%s cannot be read as a FITS file"), file);
		return 1;
	}
	int ret = 1;
	if (rx != f->rx || ry != f->ry)
		*why = g_strdup_printf(_("%s describes a %dx%d image and the grid is %dx%d"), file, rx, ry, f->rx, f->ry);
	else if (uid && strcmp(uid, f->uid))
		*why = g_strdup_printf(_("%s was exported from another Living Stack (%s)"), file, uid);
	else if (read_fits_metadata_from_path_first_HDU(file, out))
		*why = g_strdup_printf(_("the header of %s cannot be read"), file);
	else if (!has_wcs(out))
		*why = g_strdup_printf(_("%s has no WCS: plate-solve it first"), file);
	else ret = 0;
	/* Siril's plate solving flips an image that its solution shows upside
	 * down, and swaps its ROWORDER: the WCS then describes the flipped
	 * pixels, and is flipped back to describe the Stack's */
	if (!ret && out->keywords.row_order[0] && f->primary.keywords.row_order[0]
			&& strcmp(out->keywords.row_order, f->primary.keywords.row_order)) {
		out->ry = f->ry;
		flip_bottom_up_astrometry_data(out);
		siril_log_message(_("LSF: %s was flipped after it was exported (ROWORDER %s, the Stack's is %s): its WCS "
					"is flipped back to match the Stack\n"), file, out->keywords.row_order,
				f->primary.keywords.row_order);
	}
	g_free(uid);
	if (ret)
		clearfits(out);
	return ret;
}

/* the largest shift, in pixels, between the WCS and its linear part at the
 * corners and the middle of the edges of the grid: what the SIP terms
 * correct there */
double lsf_wcs_sip_shift(struct wcsprm *wcs, int rx, int ry) {
	if (!wcs || !wcs->lin.dispre)
		return 0.0;
	struct wcsprm *lin = wcs_copy_linear(wcs);
	if (!lin)
		return 0.0;
	double worst = 0.0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			double x = 0.5 * i * rx, y = 0.5 * j * ry, ra, dec, lx, ly;
			pix2wcs2(wcs, x, y, &ra, &dec);
			if (wcs2pix2(lin, ra, dec, &lx, &ly))
				continue;
			worst = fmax(worst, hypot(lx - x, ly - y));
		}
	wcsfree(lin);
	free(lin);
	return worst;
}

/* the same WCS, as far as a file stores it: the same sky position for every
 * point of a 3x3 grid, to 1e-4 arcsec, and SIP terms in both or neither */
gboolean lsf_wcs_same(struct wcsprm *a, struct wcsprm *b, int rx, int ry) {
	if (!a || !b)
		return a == b;
	if ((a->lin.dispre == NULL) != (b->lin.dispre == NULL))
		return FALSE;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			double ra1, dec1, ra2, dec2, x = 0.5 * i * rx, y = 0.5 * j * ry;
			pix2wcs2(a, x, y, &ra1, &dec1);
			pix2wcs2(b, x, y, &ra2, &dec2);
			double dra = fmod(fabs(ra1 - ra2), 360.0);
			dra = fmin(dra, 360.0 - dra) * cos(dec1 * M_PI / 180.0);
			if (hypot(dra, dec1 - dec2) * 3600.0 > 1e-4)
				return FALSE;
		}
	return TRUE;
}

/* centre, pixel scale and SIP order, for the summaries */
gchar *lsf_wcs_describe(struct wcsprm *wcs, int rx, int ry) {
	double ra, dec, cd[NAXIS][NAXIS];
	center2wcs2(wcs, rx, ry, &ra, &dec);
	wcs_cd2mat(wcs, cd);
	double scale = sqrt(fabs(cd[0][0] * cd[1][1] - cd[0][1] * cd[1][0])) * 3600.0;
	int order = wcs->lin.dispre ? extract_SIP_order_and_matrices(wcs->lin.dispre, NULL, NULL, NULL, NULL) : 0;
	if (order > 1)
		return g_strdup_printf(_("centre RA %.5f Dec %+.5f, %.3f\"/px, SIP order %d"), ra, dec, scale, order);
	return g_strdup_printf(_("centre RA %.5f Dec %+.5f, %.3f\"/px, linear"), ra, dec, scale);
}

/* The transform of a frame onto the grid from the two WCS (Section 10.4),
 * in Siril's registration convention, as calcH_from_corners (astrometric.c)
 * computes it for astrometric registration, but fitted on a grid of points
 * over the frame instead of its four corners: with SIP terms the map is not
 * exactly a homography, and the residual is spread over the frame. Between
 * two linear WCS it is exact. Returns 0, or 1 when it cannot be computed. */
int lsf_wcs_homography(struct wcsprm *img, int rx, int ry, struct wcsprm *ref, int ry_ref,
		transformation_type type, Homography *H) {
	int n = 0;
	s_star *a = calloc(LSF_WCS_GRID * LSF_WCS_GRID, sizeof(s_star));
	s_star *b = calloc(LSF_WCS_GRID * LSF_WCS_GRID, sizeof(s_star));
	if (!a || !b) {
		PRINT_ALLOC_ERR;
		free(a);
		free(b);
		return 1;
	}
	for (int i = 0; i < LSF_WCS_GRID; i++)
		for (int j = 0; j < LSF_WCS_GRID; j++) {
			/* Siril coordinates, then opencv's, as calcH_from_corners */
			double x = (double) rx * i / (LSF_WCS_GRID - 1), y = (double) ry * j / (LSF_WCS_GRID - 1);
			double ra, dec, xr, yr;
			pix2wcs2(img, x, y, &ra, &dec);
			if (wcs2pix2(ref, ra, dec, &xr, &yr))
				continue;
			a[n].x = x - 0.5;
			a[n].y = ry - y - 0.5;
			b[n].x = xr - 0.5;
			b[n].y = ry_ref - yr - 0.5;
			a[n].index = b[n].index = n;
			n++;
		}
	int ret = 1;
	if (n >= 4) {
		unsigned char *mask = cvCalculH(a, b, n, H, type, 0.f);
		ret = mask == NULL;
		free(mask);
	}
	free(a);
	free(b);
	return ret;
}

static void add_card(GString *out, const char *key, const char *value) {
	char card[FLEN_CARD];
	int status = 0;
	if (ffmkky(key, (char *) value, NULL, card, &status))
		return;
	size_t n = strlen(card);
	g_string_append_len(out, card, n);
	for (; n < 80; n++)
		g_string_append_c(out, ' ');
}

static void add_matrix(GString *out, const char *name, double m[MAX_SIP_SIZE][MAX_SIP_SIZE], int order) {
	char key[FLEN_KEYWORD], value[FLEN_VALUE];
	gboolean any = FALSE;
	for (int p = 0; p <= order; p++)
		for (int q = 0; p + q <= order; q++)
			any |= m[p][q] != 0.0;
	if (!any)
		return;
	g_snprintf(key, sizeof key, "%s_ORDER", name);
	g_snprintf(value, sizeof value, "%d", order);
	add_card(out, key, value);
	for (int p = 0; p <= order; p++)
		for (int q = 0; p + q <= order; q++) {
			if (m[p][q] == 0.0)
				continue;
			g_snprintf(key, sizeof key, "%s_%d_%d", name, p, q);
			g_snprintf(value, sizeof value, "%.15E", m[p][q]);
			add_card(out, key, value);
		}
}

/* the SIP terms of wcs and the reference pixel they apply to, as cards */
static gchar *sip_cards(const struct wcsprm *wcs) {
	double A[MAX_SIP_SIZE][MAX_SIP_SIZE] = { 0 }, B[MAX_SIP_SIZE][MAX_SIP_SIZE] = { 0 };
	double AP[MAX_SIP_SIZE][MAX_SIP_SIZE] = { 0 }, BP[MAX_SIP_SIZE][MAX_SIP_SIZE] = { 0 };
	int order = extract_SIP_order_and_matrices(wcs->lin.dispre, A, B, AP, BP);
	GString *s = g_string_new(NULL);
	char value[FLEN_VALUE];
	g_snprintf(value, sizeof value, "%.15E", wcs->crpix[0]);
	add_card(s, "CRPIX1", value);
	g_snprintf(value, sizeof value, "%.15E", wcs->crpix[1]);
	add_card(s, "CRPIX2", value);
	add_matrix(s, "A", A, order);
	add_matrix(s, "B", B, order);
	add_matrix(s, "AP", AP, order);
	add_matrix(s, "BP", BP, order);
	return g_string_free(s, FALSE);
}

/* REG_DIST of a frame (Section 9): the SIP terms of the model it was
 * undistorted with, and the reference pixel they apply to, as 80-character
 * cards. The master pattern is resolved with the frame's header, as Siril's
 * apply does; with the frames' own models (astrometric registration), frame
 * is the source frame, whose WCS holds them. NULL when the model cannot be
 * read, or when the frame has no SIP terms. */
gchar *lsf_disto_cards(const disto_params *d, fits *frame) {
	if (d && d->index == DISTO_FILES)
		return has_wcs(frame) && frame->keywords.wcslib->lin.dispre ? sip_cards(frame->keywords.wcslib) : NULL;
	if (!d || !d->filename || (d->index != DISTO_FILE && d->index != DISTO_MASTER))
		return NULL;
	int status = 0;
	gchar *path = d->index == DISTO_MASTER ? path_parse(frame, d->filename, PATHPARSE_MODE_READ, &status)
		: g_strdup(d->filename);
	if (status || !path) {
		g_free(path);
		return NULL;
	}
	fits m = { 0 };
	gchar *out = NULL;
	if (!read_fits_metadata_from_path_first_HDU(path, &m) && has_wcs(&m) && m.keywords.wcslib->lin.dispre)
		out = sip_cards(m.keywords.wcslib);
	clearfits(&m);
	g_free(path);
	return out;
}

/***** registration through the WCS (lsf_append -wcsreg) *****/

struct wcsreg_data {
	struct wcsprm *ref;
	int ref_ry, layer;
	transformation_type type;
	gboolean linear;
};

/* each frame: its transform onto the grid from its WCS; its star measures
 * are made on the grid when it is merged (Section 10.3) */
static int wcsreg_image_hook(struct generic_seq_args *args, int o, int i, fits *fit, rectangle *_, int threads) {
	struct wcsreg_data *d = args->user;
	regdata *rd = &args->seq->regparam[d->layer][i];
	int filenum = args->seq->imgparam[i].filenum;
	memset(rd, 0, sizeof(regdata));
	if (!has_wcs(fit)) {
		siril_log_message(_("LSF: image %d has not been plate-solved and is not merged\n"), filenum);
		args->seq->imgparam[i].incl = FALSE;
		return 1;
	}
	/* wcslib computes in work arrays of the wcsprm (the SIP terms in
	 * particular): frames processed in parallel each use a copy of the
	 * Stack's */
	struct wcsprm *ref = wcs_deepcopy(d->ref, NULL);
	/* on an undistorted grid the frame is undistorted by the model of the
	 * file, so only the linear part of its WCS maps it */
	struct wcsprm *w = d->linear ? wcs_copy_linear(fit->keywords.wcslib) : fit->keywords.wcslib;
	int ret = !ref || !w || lsf_wcs_homography(w, fit->rx, fit->ry, ref, d->ref_ry, d->type, &rd->H);
	if (d->linear && w) {
		wcsfree(w);
		free(w);
	}
	if (ref) {
		wcsfree(ref);
		free(ref);
	}
	if (ret) {
		siril_log_message(_("LSF: the WCS of image %d cannot be mapped onto the Stack, it is not merged\n"), filenum);
		SetNullH(&rd->H);
		args->seq->imgparam[i].incl = FALSE;
		return 1;
	}
	return 0;
}

/* Fills seq->regparam[layer] with the transforms onto the grid of the
 * selected frames, computed from their WCS and the Stack's (ref, of height
 * ref_ry). Frames that are not plate-solved are unselected. linear: map
 * only the linear part of the frames' WCS, on an undistorted grid. Returns
 * 0, or 1 when no frame could be registered. */
int lsf_register_by_wcs(sequence *seq, int layer, struct wcsprm *ref, int ref_ry, transformation_type type,
		gboolean linear, int max_threads) {
	if (!seq->regparam[layer]) {
		seq->regparam[layer] = calloc(seq->number, sizeof(regdata));
		if (!seq->regparam[layer]) {
			PRINT_ALLOC_ERR;
			return 1;
		}
	}
	struct wcsreg_data d = { 0 };
	d.ref = ref;
	d.ref_ry = ref_ry;
	d.layer = layer;
	d.type = type;
	d.linear = linear;

	struct generic_seq_args *args = create_default_seqargs(seq);
	args->filtering_criterion = seq_filter_included;
	args->nb_filtered_images = seq->selnum;
	args->image_hook = wcsreg_image_hook;
	args->stop_on_error = FALSE;
	args->description = _("LSF registration through the WCS");
	args->has_output = FALSE;
	args->already_in_a_thread = TRUE;
	args->max_parallel_images = MAX(max_threads, 1);
	args->user = &d;
	generic_sequence_worker(args);
	free(args);
	fix_selnum(seq, FALSE);
	int registered = 0;
	for (int i = 0; i < seq->number; i++)
		if (seq->imgparam[i].incl && guess_transform_from_H(seq->regparam[layer][i].H) != NULL_TRANSFORMATION)
			registered++;
	return registered == 0;
}
