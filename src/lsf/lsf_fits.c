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

/* LSF file I/O (FITS profile v0.4): validation (Section 8), writing
 * (Sections 3 to 9) and the commit protocol (Section 10.2, 10.11) */

#include <string.h>
#include <limits.h>
#include <math.h>
#ifdef _WIN32
/* must come before fitsio.h: it #defines TBYTE, which clashes with winnt.h */
#include <windows.h>
#include <io.h>
#else
#include <unistd.h>
#endif
#include <fitsio.h>
#include <glib/gstdio.h>
#include <fcntl.h>

#include "core/siril.h"
#include "core/siril_log.h"
#include "io/image_format_fits.h"
#include "lsf/lsf.h"

static const char *plane_kind[LSF_NB_PL] = { "MEAN", "W", "M2", "W2" };
static const char *chan_suffix[3] = { "_R", "_G", "_B" };

/* keywords that only make sense in an LSF file, removed by lsf_export */
static const char *lsf_keys[] = { "LSF", "LSFVER", "LSFUID", "LSFGEN", "LSFSENS", "LSFBAYR",
	"LSFDMOD", "LSFNORM", "LSFNMIN", "LSFNMAX", "LSFHASH", "LSFBUSY" };

/* must never appear in the primary (Section 4.1, 7) */
static const char *primary_banned_keys[] = { "BAYERPAT", "XBAYROFF", "YBAYROFF",
	"MIPS-HI", "MIPS-LO", "MIPS-FHI", "MIPS-FLO" };

int lsf_history_alloc(lsf_history *h, int rx, int ry, int nchan) {
	memset(h, 0, sizeof(lsf_history));
	h->rx = rx; h->ry = ry; h->nchan = nchan;
	size_t n = (size_t) rx * ry;
	for (int c = 0; c < nchan; c++) {
		for (int k = 0; k < LSF_NB_PL; k++) {
			h->plane[k][c] = calloc(n, sizeof(float));
			if (!h->plane[k][c]) goto fail;
		}
		h->count[c] = calloc(n, sizeof(WORD));
		if (!h->count[c]) goto fail;
	}
	return 0;
fail:
	PRINT_ALLOC_ERR;
	lsf_history_free(h);
	return 1;
}

void lsf_history_free(lsf_history *h) {
	for (int c = 0; c < 3; c++) {
		for (int k = 0; k < LSF_NB_PL; k++) {
			free(h->plane[k][c]);
			h->plane[k][c] = NULL;
		}
		free(h->count[c]);
		h->count[c] = NULL;
	}
}

void lsf_file_free(lsf_file *f) {
	g_free(f->path);
	f->path = NULL;
	lsf_settings_free(&f->settings);
	lsf_history_free(&f->hist);
	clearfits(&f->primary);
	free(f->stack_data);
	f->stack_data = NULL;
	if (f->unknown_hdus)
		g_array_free(f->unknown_hdus, TRUE);
	f->unknown_hdus = NULL;
}

/* SHA-256 of float32 data as stored in FITS: big-endian, file order, no padding */
gchar *lsf_hash_float_data(const float *data, size_t n) {
	GChecksum *ck = g_checksum_new(G_CHECKSUM_SHA256);
	guint32 buf[4096];
	size_t done = 0;
	while (done < n) {
		size_t chunk = MIN(n - done, G_N_ELEMENTS(buf));
		for (size_t i = 0; i < chunk; i++) {
			guint32 v;
			memcpy(&v, data + done + i, sizeof v);
			buf[i] = GUINT32_TO_BE(v);
		}
		g_checksum_update(ck, (const guchar *) buf, chunk * sizeof(guint32));
		done += chunk;
	}
	gchar *hex = g_strdup(g_checksum_get_string(ck));	// lowercase hex
	g_checksum_free(ck);
	return hex;
}

static gchar *plane_name(int kind, int chan, int nchan) {
	if (nchan == 1)
		return g_strdup(plane_kind[kind]);
	return g_strdup_printf("%s%s", plane_kind[kind], chan_suffix[chan]);
}

/***************************** reading *****************************/

static int validation_failure(int check, const char *path, const char *what) {
	siril_log_error(_("LSF validation of %s failed at check %d (Section 8): %s. "
				"The file was not modified.\n"), path, check, what);
	return check;
}

/* reads LSFUID and LSFGEN of the current HDU, returns non-zero if missing */
static int read_uid_gen(fitsfile *fptr, char *uid, int *gen) {
	int status = 0;
	char value[FLEN_VALUE];
	fits_read_key(fptr, TSTRING, "LSFUID", value, NULL, &status);
	fits_read_key(fptr, TINT, "LSFGEN", gen, NULL, &status);
	if (!status)
		g_strlcpy(uid, value, FLEN_VALUE);
	return status;
}

static GPtrArray *read_settings_rows(fitsfile *fptr, int *status) {
	long nrows = 0;
	int ckey, ctype, cvalue;
	fits_get_num_rows(fptr, &nrows, status);
	fits_get_colnum(fptr, CASEINSEN, "KEY", &ckey, status);
	fits_get_colnum(fptr, CASEINSEN, "TYPE", &ctype, status);
	fits_get_colnum(fptr, CASEINSEN, "VALUE", &cvalue, status);
	if (*status)
		return NULL;
	GPtrArray *rows = g_ptr_array_new_with_free_func(lsf_setting_row_free);
	for (long r = 1; r <= nrows && !*status; r++) {
		char key[80] = { 0 }, type[80] = { 0 };
		char *pk = key, *pt = type;
		long repeat = 0, offset = 0;
		fits_read_col_str(fptr, ckey, r, 1, 1, "", &pk, NULL, status);
		fits_read_col_str(fptr, ctype, r, 1, 1, "", &pt, NULL, status);
		fits_read_descript(fptr, cvalue, r, &repeat, &offset, status);
		char *value = calloc(repeat + 1, 1);
		fits_read_col_str(fptr, cvalue, r, 1, 1, "", &value, NULL, status);
		lsf_setting_row *row = g_new(lsf_setting_row, 1);
		row->key = g_strstrip(g_strdup(key));
		row->type = g_strstrip(g_strdup(type));
		row->value = g_strdup(value);
		free(value);
		g_ptr_array_add(rows, row);
	}
	return rows;
}

static const char *find_row(GPtrArray *rows, const char *key) {
	for (guint i = 0; i < rows->len; i++) {
		lsf_setting_row *r = g_ptr_array_index(rows, i);
		if (!strcmp(r->key, key))
			return r->value;
	}
	return NULL;
}

/* Opens an LSF file, runs the pre-merge validation of Section 8 in order,
 * and loads what is needed to append to it. Returns 0 on success, or the
 * number of the failed check. The file is never modified.
 * With tolerate_edited_stack, a failure of check 8 alone (the Stack was
 * modified outside Siril, the accumulator is intact) is not an error: the
 * file is loaded, f->stack_modified is set and f->stack_data holds the
 * edited Stack, for Repair and Adopt. */
int lsf_open_and_validate(const char *path, lsf_file *f, gboolean tolerate_edited_stack) {
	fitsfile *fptr = NULL;
	int status = 0, hdunum = 0;
	int ret = 0;
	memset(f, 0, sizeof(lsf_file));
	f->path = g_strdup(path);
	f->unknown_hdus = g_array_new(FALSE, FALSE, sizeof(int));
	lsf_settings_init(&f->settings);
	GPtrArray *rows = NULL;
	float *primary_data = NULL;

	if (fits_open_diskfile(&fptr, path, READONLY, &status)) {
		report_fits_error(status);
		return validation_failure(1, path, _("cannot open the file"));
	}

	/* check 1: LSF = T and supported version */
	int lsf = 0;
	char ver[FLEN_VALUE] = { 0 };
	status = 0;
	fits_read_key(fptr, TLOGICAL, "LSF", &lsf, NULL, &status);
	fits_read_key(fptr, TSTRING, "LSFVER", ver, NULL, &status);
	if (status || !lsf) {
		ret = validation_failure(1, path, _("not an LSF file (LSF keyword missing or false)"));
		goto end;
	}
	if (strcmp(ver, LSF_VERSION)) {
		if (!strcmp(ver, "0.3"))
			ret = validation_failure(1, path, _("LSF version 0.3 files must be converted to 0.4 first"));
		else ret = validation_failure(1, path, _("unsupported LSF version"));
		goto end;
	}

	/* check 2: not busy */
	int busy = 0;
	status = 0;
	fits_read_key(fptr, TLOGICAL, "LSFBUSY", &busy, NULL, &status);
	if (!status && busy) {
		ret = validation_failure(2, path, _("LSFBUSY is set: an update was interrupted, restore a backup"));
		goto end;
	}

	/* check 3: required HDUs, found by name, with consistent dimensions */
	status = 0;
	fits_get_num_hdus(fptr, &hdunum, &status);
	int hdu_acc = 0, hdu_count = 0, hdu_settings = 0, hdu_ledger = 0, hdu_rej = 0;
	for (int i = 2; i <= hdunum; i++) {
		char extname[FLEN_VALUE] = { 0 };
		status = 0;
		fits_movabs_hdu(fptr, i, NULL, &status);
		fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status);
		if (status) extname[0] = '\0';
		if (!strcmp(extname, LSF_EXT_ACC)) hdu_acc = i;
		else if (!strcmp(extname, LSF_EXT_COUNT)) hdu_count = i;
		else if (!strcmp(extname, LSF_EXT_SETTINGS)) hdu_settings = i;
		else if (!strcmp(extname, LSF_EXT_LEDGER)) hdu_ledger = i;
		else {
			if (!strcmp(extname, "LSF_REJ")) hdu_rej = i;
			if (!strcmp(extname, "LSF_REJ") || !strcmp(extname, "LSF_STARS")) {
				/* optional LSF HDUs from a later writer: dropped, they would be
				 * stale after this append */
				siril_log_warning(_("LSF: %s is not supported by this version and will be dropped\n"), extname);
				continue;
			}
			g_array_append_val(f->unknown_hdus, i);
		}
	}
	if (!hdu_acc || !hdu_count || !hdu_settings || !hdu_ledger) {
		ret = validation_failure(3, path, _("a required LSF HDU is missing"));
		goto end;
	}
	f->ledger_hdu = hdu_ledger;

	int naxis = 0, bitpix = 0;
	long naxes[3] = { 0, 0, 1 };
	status = 0;
	fits_movabs_hdu(fptr, 1, NULL, &status);
	fits_get_img_param(fptr, 3, &bitpix, &naxis, naxes, &status);
	if (status || bitpix != FLOAT_IMG || (naxis != 2 && !(naxis == 3 && naxes[2] == 3))) {
		ret = validation_failure(3, path, _("the primary is not a float32 mono or RGB image"));
		goto end;
	}
	f->rx = naxes[0];
	f->ry = naxes[1];
	f->nchan = naxis == 3 ? 3 : 1;

	long anaxes[3] = { 0 };
	int anaxis = 0, abitpix = 0;
	status = 0;
	fits_movabs_hdu(fptr, hdu_acc, NULL, &status);
	fits_get_img_param(fptr, 3, &abitpix, &anaxis, anaxes, &status);
	if (status || abitpix != FLOAT_IMG || anaxis != 3 || anaxes[0] != f->rx || anaxes[1] != f->ry
			|| anaxes[2] != LSF_NB_PL * f->nchan) {
		ret = validation_failure(3, path, _("LSF_ACC does not have the expected dimensions"));
		goto end;
	}
	/* map the planes by name (C6) */
	int plane_index[LSF_NB_PL][3];
	for (int k = 0; k < LSF_NB_PL; k++)
		for (int c = 0; c < f->nchan; c++) {
			plane_index[k][c] = -1;
			gchar *expected = plane_name(k, c, f->nchan);
			for (int p = 0; p < anaxes[2]; p++) {
				char key[FLEN_KEYWORD], val[FLEN_VALUE] = { 0 };
				g_snprintf(key, sizeof key, "LSFPL%d", p + 1);
				status = 0;
				fits_read_key(fptr, TSTRING, key, val, NULL, &status);
				if (!status && !strcmp(g_strstrip(val), expected))
					plane_index[k][c] = p;
			}
			g_free(expected);
			if (plane_index[k][c] < 0) {
				if (k == LSF_PL_W2)
					ret = validation_failure(3, path, _("the file has no W2 plane, appending is refused (Section 11)"));
				else ret = validation_failure(3, path, _("LSF_ACC does not hold the expected named planes"));
				goto end;
			}
		}

	long cnaxes[3] = { 0, 0, 1 };
	int cnaxis = 0, cbitpix = 0;
	status = 0;
	fits_movabs_hdu(fptr, hdu_count, NULL, &status);
	fits_get_img_param(fptr, 3, &cbitpix, &cnaxis, cnaxes, &status);
	if (status || cnaxes[0] != f->rx || cnaxes[1] != f->ry ||
			(f->nchan == 1 && cnaxis != 2) ||
			(f->nchan == 3 && !(cnaxis == 3 && (cnaxes[2] == 3 || cnaxes[2] == 1)))) {
		ret = validation_failure(3, path, _("LSF_COUNT does not have the expected dimensions"));
		goto end;
	}

	/* check 4: LSFUID and LSFGEN identical across all LSF HDUs */
	int hdus[] = { 1, hdu_acc, hdu_count, hdu_settings, hdu_ledger };
	char uid0[FLEN_VALUE] = { 0 };
	int gen0 = -1;
	for (int i = 0; i < (int) G_N_ELEMENTS(hdus); i++) {
		char uid[FLEN_VALUE] = { 0 };
		int gen = -1;
		status = 0;
		fits_movabs_hdu(fptr, hdus[i], NULL, &status);
		if (status || read_uid_gen(fptr, uid, &gen)) {
			ret = validation_failure(4, path, _("LSFUID or LSFGEN missing in an LSF HDU"));
			goto end;
		}
		if (i == 0) {
			g_strlcpy(uid0, uid, sizeof uid0);
			gen0 = gen;
		} else if (strcmp(uid, uid0) || gen != gen0) {
			ret = validation_failure(4, path, _("LSFUID or LSFGEN differ between LSF HDUs"));
			goto end;
		}
	}
	g_strlcpy(f->uid, uid0, sizeof f->uid);
	f->gen = gen0;

	/* check 5: DATASUM and CHECKSUM on every extension where present */
	for (int i = 2; i <= hdunum; i++) {
		int dataok = 0, hduok = 0;
		status = 0;
		fits_movabs_hdu(fptr, i, NULL, &status);
		fits_verify_chksum(fptr, &dataok, &hduok, &status);
		if (status || dataok == -1 || hduok == -1) {
			ret = validation_failure(5, path, _("a checksum does not match, the accumulator cannot be trusted"));
			goto end;
		}
	}

	/* check 6: ledger generations */
	status = 0;
	fits_movabs_hdu(fptr, hdu_ledger, NULL, &status);
	fits_get_num_rows(fptr, &f->ledger_rows, &status);
	int cgen = 0, cbatch = 0;
	fits_get_colnum(fptr, CASEINSEN, "GEN", &cgen, &status);
	if (status) {
		ret = validation_failure(6, path, _("the ledger has no GEN column"));
		goto end;
	}
	f->max_ledger_gen = 0;
	f->nb_batches = 0;
	int has_batch = !fits_get_colnum(fptr, CASEINSEN, "BATCH", &cbatch, &status);
	status = 0;
	for (long r = 1; r <= f->ledger_rows; r++) {
		int g = 0, b = 0;
		fits_read_col(fptr, TINT, cgen, r, 1, 1, NULL, &g, NULL, &status);
		if (has_batch)
			fits_read_col(fptr, TINT, cbatch, r, 1, 1, NULL, &b, NULL, &status);
		if (status) break;
		if (g > f->max_ledger_gen) f->max_ledger_gen = g;
		if (b > f->nb_batches) f->nb_batches = b;
	}
	/* a generation without ledger rows changed no frame (Repair, Adopt) */
	if (status || f->max_ledger_gen > f->gen) {
		ret = validation_failure(6, path, _("ledger generations do not match LSFGEN"));
		goto end;
	}

	/* check 7: mirrored settings agree */
	status = 0;
	fits_movabs_hdu(fptr, hdu_settings, NULL, &status);
	rows = read_settings_rows(fptr, &status);
	if (status || !rows) {
		ret = validation_failure(7, path, _("LSF_SETTINGS cannot be read"));
		goto end;
	}
	if (lsf_settings_from_rows(&f->settings, rows, f->nchan)) {
		ret = validation_failure(7, path, _("LSF_SETTINGS holds unsupported or malformed values"));
		goto end;
	}
	{
		char sens[FLEN_VALUE] = { 0 }, bayr[FLEN_VALUE] = { 0 }, dmod[FLEN_VALUE] = { 0 };
		int norm = 0;
		double nmin = 0.0, nmax = 0.0;
		int nstatus = 0;
		status = 0;
		fits_movabs_hdu(fptr, 1, NULL, &status);
		fits_read_key(fptr, TSTRING, "LSFSENS", sens, NULL, &status);
		fits_read_key(fptr, TSTRING, "LSFDMOD", dmod, NULL, &status);
		fits_read_key(fptr, TLOGICAL, "LSFNORM", &norm, NULL, &status);
		const char *s_class = find_row(rows, "sensor.class");
		const char *s_dmod = find_row(rows, "cfa.mode");
		gboolean bad = status || !s_class || !s_dmod || g_ascii_strcasecmp(sens, s_class)
			|| g_ascii_strcasecmp(dmod, s_dmod) || (!norm) != (!f->settings.outnorm);
		if (!bad && f->settings.is_cfa) {
			nstatus = 0;
			fits_read_key(fptr, TSTRING, "LSFBAYR", bayr, NULL, &nstatus);
			bad = nstatus || g_ascii_strcasecmp(bayr, f->settings.bayer);
		}
		if (!bad && norm) {
			nstatus = 0;
			fits_read_key(fptr, TDOUBLE, "LSFNMIN", &nmin, NULL, &nstatus);
			fits_read_key(fptr, TDOUBLE, "LSFNMAX", &nmax, NULL, &nstatus);
			/* Section 4.2: refuse to work on a normalized primary without them */
			bad = nstatus || nmin != f->settings.outnorm_min || nmax != f->settings.outnorm_max;
			f->outnorm_min = nmin;
			f->outnorm_max = nmax;
		}
		if (!bad && (f->settings.rejmaps != (hdu_rej != 0)))
			bad = TRUE;
		if (bad) {
			ret = validation_failure(7, path, _("primary header keywords disagree with LSF_SETTINGS"));
			goto end;
		}
	}

	/* check 8: the primary data still matches LSFHASH */
	{
		char hash[FLEN_VALUE] = { 0 };
		size_t n = (size_t) f->rx * f->ry * f->nchan;
		primary_data = malloc(n * sizeof(float));
		if (!primary_data) {
			PRINT_ALLOC_ERR;
			ret = 8;
			goto end;
		}
		long fpixel[3] = { 1, 1, 1 };
		status = 0;
		fits_movabs_hdu(fptr, 1, NULL, &status);
		fits_read_key(fptr, TSTRING, "LSFHASH", hash, NULL, &status);
		fits_read_pix(fptr, TFLOAT, fpixel, n, NULL, primary_data, NULL, &status);
		gchar *computed = status ? NULL : lsf_hash_float_data(primary_data, n);
		gboolean match = computed && !g_ascii_strcasecmp(computed, hash);
		g_free(computed);
		if (!match) {
			if (!tolerate_edited_stack) {
				ret = validation_failure(8, path, _("the Stack (HDU 0) was modified outside Siril. "
							"Run lsf_check with -repair to rebuild it from the accumulator, "
							"or with -adopt to keep the edit"));
				goto end;
			}
			siril_log_warning(_("LSF validation of %s: check 8 (Section 8) failed, the Stack (HDU 0) "
						"was modified outside Siril. The accumulator is intact.\n"), path);
			f->stack_modified = TRUE;
			f->stack_data = primary_data;
			primary_data = NULL;
		}
	}

	/* load the accumulator */
	if (lsf_history_alloc(&f->hist, f->rx, f->ry, f->nchan)) {
		ret = -1;
		goto end;
	}
	status = 0;
	fits_movabs_hdu(fptr, hdu_acc, NULL, &status);
	size_t npix = (size_t) f->rx * f->ry;
	for (int k = 0; k < LSF_NB_PL && !status; k++)
		for (int c = 0; c < f->nchan && !status; c++) {
			long fpixel[3] = { 1, 1, plane_index[k][c] + 1 };
			fits_read_pix(fptr, TFLOAT, fpixel, npix, NULL, f->hist.plane[k][c], NULL, &status);
		}
	fits_movabs_hdu(fptr, hdu_count, NULL, &status);
	for (int c = 0; c < f->nchan && !status; c++) {
		/* a single count plane is shared by all channels (Section 3) */
		long fpixel[3] = { 1, 1, (cnaxis == 3 && cnaxes[2] == 3) ? c + 1 : 1 };
		fits_read_pix(fptr, TUSHORT, fpixel, npix, NULL, f->hist.count[c], NULL, &status);
	}
	if (status) {
		report_fits_error(status);
		ret = -1;
		goto end;
	}
	{
		int sat = 0, st = 0;
		fits_read_key(fptr, TLOGICAL, "LSFCSAT", &sat, NULL, &st);
		f->hist.count_saturated = !st && sat;
	}

end:
	free(primary_data);
	if (rows)
		g_ptr_array_free(rows, TRUE);
	status = 0;
	if (fptr)
		fits_close_file(fptr, &status);
	if (!ret) {
		/* header metadata of the primary, for the next write */
		if (read_fits_metadata_from_path(path, &f->primary)) {
			siril_log_error(_("LSF: could not read the primary header of %s\n"), path);
			ret = -1;
		}
	}
	return ret;
}

/***************************** writing *****************************/

static gboolean is_wcs_key(const char *key) {
	static const char *prefixes[] = { "CTYPE", "CRVAL", "CRPIX", "CDELT", "CROTA", "CUNIT",
		"CD1_", "CD2_", "PC1_", "PC2_", "PV1_", "PV2_", "A_", "B_", "AP_", "BP_", "WCSAXES",
		"LONPOLE", "LATPOLE", "RADESYS", "EQUINOX", "MJDREF", "WCSNAME" };
	for (int i = 0; i < (int) G_N_ELEMENTS(prefixes); i++)
		if (g_str_has_prefix(key, prefixes[i]))
			return TRUE;
	return FALSE;
}

/* deletes a keyword that may be absent, without leaving cfitsio error
 * messages that a later report_fits_error() would print */
static void delete_key_quiet(fitsfile *fptr, const char *key) {
	int st = 0;
	fits_write_errmark();
	fits_delete_key(fptr, key, &st);
	fits_clear_errmark();
}

/* the WCS cards of the current HDU, as records */
static GPtrArray *collect_wcs_cards(fitsfile *fptr) {
	GPtrArray *cards = g_ptr_array_new_with_free_func(g_free);
	int nkeys = 0, status = 0;
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);
	for (int i = 1; i <= nkeys && !status; i++) {
		char card[FLEN_CARD], key[FLEN_KEYWORD];
		int len = 0;
		fits_read_record(fptr, i, card, &status);
		if (status) break;
		fits_get_keyname(card, key, &len, &status);
		/* for a cube the WCS describes axes 1 and 2 only (Section 5.1) */
		gboolean third_axis = strlen(key) > 1 && key[strlen(key) - 1] == '3' && key[0] == 'C';
		if (!status && is_wcs_key(key) && !third_axis && strcmp(key, "WCSAXES"))
			g_ptr_array_add(cards, g_strdup(card));
		status = 0;
	}
	return cards;
}

static void write_cards(fitsfile *fptr, GPtrArray *cards, int *status) {
	for (guint i = 0; i < cards->len; i++)
		fits_write_record(fptr, g_ptr_array_index(cards, i), status);
}

static void write_identity(fitsfile *fptr, const char *extname, const char *uid, int gen, int *status) {
	if (extname) {
		fits_update_key(fptr, TSTRING, "EXTNAME", (void *) extname, NULL, status);
		int one = 1;
		fits_update_key(fptr, TINT, "EXTVER", &one, NULL, status);
	}
	fits_update_key(fptr, TSTRING, "LSFUID", (void *) uid, "LSF unique id", status);
	fits_update_key(fptr, TINT, "LSFGEN", &gen, "LSF generation", status);
}

static int write_settings_hdu(fitsfile *fptr, const lsf_settings *s, int nchan,
		const char *uid, int gen, int *status) {
	char *ttype[] = { "KEY", "TYPE", "VALUE" };
	char *tform[] = { "32A", "8A", "1PA" };
	GPtrArray *rows = lsf_settings_to_rows(s, nchan);
	fits_create_tbl(fptr, BINARY_TBL, 0, 3, ttype, tform, NULL, LSF_EXT_SETTINGS, status);
	for (guint i = 0; i < rows->len && !*status; i++) {
		lsf_setting_row *r = g_ptr_array_index(rows, i);
		char *k = r->key, *t = r->type, *v = r->value;
		fits_write_col_str(fptr, 1, i + 1, 1, 1, &k, status);
		fits_write_col_str(fptr, 2, i + 1, 1, 1, &t, status);
		fits_write_col_str(fptr, 3, i + 1, 1, 1, &v, status);
	}
	write_identity(fptr, NULL, uid, gen, status);
	g_ptr_array_free(rows, TRUE);
	return *status;
}

/* ledger columns, Section 9. A '%d' in the form is the plane count P. */
static const struct { const char *name, *form; } ledger_cols[] = {
	{ "FRAME_ID", "32A" }, { "BATCH", "J" }, { "GEN", "J" }, { "ADDED_UTC", "24A" },
	{ "FILENAME", "1PA" }, { "DATE_OBS", "24A" }, { "EXPTIME", "D" }, { "GAIN", "J" },
	{ "OFFSET", "J" }, { "XBINNING", "J" }, { "YBINNING", "J" }, { "CCD_TEMP", "E" },
	{ "FILTER", "16A" }, { "ORIGHDR", "1PA" }, { "REG_TYPE", "8A" }, { "REG_XFRM", "9D" },
	{ "ALIGN", "1PA" }, { "NORM_SC", "%dD" }, { "NORM_OFF", "%dD" }, { "NORM_MTH", "%dA" },
	{ "WEIGHT", "%dD" }, { "FWHM", "E" }, { "WFWHM", "E" }, { "ROUND", "E" }, { "REGQ", "E" },
	{ "REGBG", "E" }, { "NSTARS", "J" }, { "BKG_MED", "%dE" }, { "NOISE", "%dE" },
	{ "NOISE_N", "%dE" }, { "STAT_MEAN", "%dE" }, { "STAT_SIG", "%dE" }, { "STAT_MAD", "%dE" },
	{ "STAT_SQBW", "%dE" }, { "STAT_LOC", "%dE" }, { "STAT_SCL", "%dE" }
};

static int create_ledger_hdu(fitsfile *fptr, int nchan, int *status) {
	int n = G_N_ELEMENTS(ledger_cols);
	char **ttype = g_new(char *, n), **tform = g_new(char *, n);
	for (int i = 0; i < n; i++) {
		ttype[i] = (char *) ledger_cols[i].name;
		if (!strcmp(ledger_cols[i].name, "NORM_MTH"))
			tform[i] = g_strdup_printf("%dA", 12 * nchan);
		else if (strchr(ledger_cols[i].form, '%'))
			tform[i] = g_strdup_printf(ledger_cols[i].form, nchan);
		else tform[i] = g_strdup(ledger_cols[i].form);
	}
	fits_create_tbl(fptr, BINARY_TBL, 0, n, ttype, tform, NULL, LSF_EXT_LEDGER, status);
	for (int i = 0; i < n; i++) {
		int col = i + 1;
		char key[FLEN_KEYWORD];
		if (!strcmp(tform[i], "J")) {
			long tnull = INT_MIN;
			g_snprintf(key, sizeof key, "TNULL%d", col);
			fits_update_key(fptr, TLONG, key, &tnull, NULL, status);
		}
		if (!strcmp(ledger_cols[i].name, "REG_XFRM")) {
			long dims[2] = { 3, 3 };
			fits_write_tdim(fptr, col, 2, dims, status);
		}
		if (!strcmp(ledger_cols[i].name, "NORM_MTH")) {
			long dims[2] = { 12, nchan };
			fits_write_tdim(fptr, col, 2, dims, status);
		}
		g_free(tform[i]);
	}
	g_free(ttype);
	g_free(tform);
	return *status;
}

static int colnum(fitsfile *fptr, const char *name) {
	int status = 0, c = 0;
	fits_get_colnum(fptr, CASEINSEN, (char *) name, &c, &status);
	return status ? 0 : c;
}

static void put_str(fitsfile *fptr, const char *name, long row, const char *v, int *status) {
	int c = colnum(fptr, name);
	if (!c) return;
	char *s = (char *) (v ? v : "");
	fits_write_col_str(fptr, c, row, 1, 1, &s, status);
}

static void put(fitsfile *fptr, int type, const char *name, long row, long n, void *v, int *status) {
	int c = colnum(fptr, name);
	if (!c) return;	// column absent in an older ledger
	fits_write_col(fptr, type, c, row, 1, n, v, status);
}

static void write_ledger_rows(fitsfile *fptr, GPtrArray *rows, long first_row, int nchan,
		int batch, int gen, const char *added, int *status) {
	for (guint i = 0; i < rows->len && !*status; i++) {
		lsf_ledger_row *r = g_ptr_array_index(rows, i);
		long row = first_row + i;
		put_str(fptr, "FRAME_ID", row, r->frame_id, status);
		put(fptr, TINT, "BATCH", row, 1, &batch, status);
		put(fptr, TINT, "GEN", row, 1, &gen, status);
		put_str(fptr, "ADDED_UTC", row, added, status);
		put_str(fptr, "FILENAME", row, r->filename, status);
		put_str(fptr, "DATE_OBS", row, r->date_obs, status);
		put(fptr, TDOUBLE, "EXPTIME", row, 1, &r->exptime, status);
		put(fptr, TINT, "GAIN", row, 1, &r->gain, status);
		put(fptr, TINT, "OFFSET", row, 1, &r->offset, status);
		put(fptr, TINT, "XBINNING", row, 1, &r->xbin, status);
		put(fptr, TINT, "YBINNING", row, 1, &r->ybin, status);
		put(fptr, TFLOAT, "CCD_TEMP", row, 1, &r->ccd_temp, status);
		put_str(fptr, "FILTER", row, r->filter, status);
		put_str(fptr, "ORIGHDR", row, r->orighdr, status);
		put_str(fptr, "REG_TYPE", row, r->reg_type, status);
		put(fptr, TDOUBLE, "REG_XFRM", row, 9, r->xfrm, status);
		put_str(fptr, "ALIGN", row, r->align_json, status);
		put(fptr, TDOUBLE, "NORM_SC", row, nchan, r->norm_sc, status);
		put(fptr, TDOUBLE, "NORM_OFF", row, nchan, r->norm_off, status);
		{
			char mth[3 * 12 + 1];
			for (int c = 0; c < nchan; c++)
				g_snprintf(mth + 12 * c, 13, "%-12s", r->norm_mth[c]);
			put_str(fptr, "NORM_MTH", row, mth, status);
		}
		put(fptr, TDOUBLE, "WEIGHT", row, nchan, r->weight, status);
		put(fptr, TFLOAT, "FWHM", row, 1, &r->fwhm, status);
		put(fptr, TFLOAT, "WFWHM", row, 1, &r->wfwhm, status);
		put(fptr, TFLOAT, "ROUND", row, 1, &r->round, status);
		put(fptr, TFLOAT, "REGQ", row, 1, &r->regq, status);
		put(fptr, TFLOAT, "REGBG", row, 1, &r->regbg, status);
		put(fptr, TINT, "NSTARS", row, 1, &r->nstars, status);
		put(fptr, TFLOAT, "BKG_MED", row, nchan, r->bkg_med, status);
		put(fptr, TFLOAT, "NOISE", row, nchan, r->noise, status);
		put(fptr, TFLOAT, "NOISE_N", row, nchan, r->noise_n, status);
		put(fptr, TFLOAT, "STAT_MEAN", row, nchan, r->st_mean, status);
		put(fptr, TFLOAT, "STAT_SIG", row, nchan, r->st_sig, status);
		put(fptr, TFLOAT, "STAT_MAD", row, nchan, r->st_mad, status);
		put(fptr, TFLOAT, "STAT_SQBW", row, nchan, r->st_sqbw, status);
		put(fptr, TFLOAT, "STAT_LOC", row, nchan, r->st_loc, status);
		put(fptr, TFLOAT, "STAT_SCL", row, nchan, r->st_scl, status);
	}
}

static int flush_to_disk(const char *path) {
	int fd = g_open(path, O_RDWR, 0);
	if (fd < 0)
		return 1;
#ifdef _WIN32
	int ret = _commit(fd);
#else
	int ret = fsync(fd);
#endif
	close(fd);
	return ret;
}

/* atomic replace of target by source (Section 10.11). On Windows a file
 * that was just written can be held for a moment by another process, such as
 * a virus scanner: a sharing or access error is retried for about 2 s before
 * giving up. */
static int replace_file(const char *source, const char *target) {
#ifdef _WIN32
	wchar_t *wsrc = g_utf8_to_utf16(source, -1, NULL, NULL, NULL);
	wchar_t *wdst = g_utf8_to_utf16(target, -1, NULL, NULL, NULL);
	DWORD err = 0;
	BOOL ok = FALSE;
	for (int attempt = 0; attempt < 20 && !ok; attempt++) {
		ok = MoveFileExW(wsrc, wdst, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH);
		if (ok) break;
		err = GetLastError();
		if (err != ERROR_SHARING_VIOLATION && err != ERROR_ACCESS_DENIED && err != ERROR_LOCK_VIOLATION)
			break;
		g_usleep(100000);
	}
	g_free(wsrc);
	g_free(wdst);
	if (!ok)
		siril_log_warning(_("LSF: replacing %s failed with Windows error %lu\n"), target, (unsigned long) err);
	return ok ? 0 : 1;
#else
	return g_rename(source, target);
#endif
}

static gchar *now_utc_iso(void) {
	GDateTime *now = g_date_time_new_now_utc();
	gchar *s = g_date_time_format(now, "%Y-%m-%dT%H:%M:%S");
	g_date_time_unref(now);
	return s;
}

/* Writes the complete new file to a temporary name, then atomically replaces
 * path with it (Section 10.2 steps 6 to 8). primary holds the header
 * metadata and the float pixels of HDU 0. old is the validated previous
 * file, or NULL on initialization. */
int lsf_write_file(const char *path, lsf_file *old, const lsf_settings *s, lsf_history *h,
		fits *primary, const char *uid, int gen, GPtrArray *new_rows, int batch, gboolean backup) {
	int status = 0;
	fitsfile *fptr = NULL, *oldf = NULL;
	gchar *tmp = g_strdup_printf("%s.lsf_tmp", path);
	gchar *added = now_utc_iso();
	GPtrArray *wcs = NULL;
	int nchan = h->nchan;
	size_t npix = (size_t) h->rx * h->ry;
	int retval = 1;

	g_unlink(tmp);
	if (siril_fits_create_diskfile(&fptr, tmp, &status)) {
		report_fits_error(status);
		goto end;
	}

	/* HDU 0: the Stack (Section 4) */
	long naxes[3] = { h->rx, h->ry, nchan };
	int naxis = nchan == 3 ? 3 : 2;
	fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
	if (status) {
		report_fits_error(status);
		goto end;
	}
	primary->fptr = fptr;
	primary->bitpix = FLOAT_IMG;
	primary->orig_bitpix = FLOAT_IMG;
	if (save_opened_fits(primary)) {
		siril_log_error(_("LSF: could not write the primary HDU\n"));
		goto end;
	}
	status = 0;
	for (int i = 0; i < (int) G_N_ELEMENTS(primary_banned_keys); i++) {
		delete_key_quiet(fptr, primary_banned_keys[i]);
	}
	{
		int t = TRUE, f = FALSE;
		double one = 1.0, zero = 0.0;
		fits_update_key(fptr, TDOUBLE, "BSCALE", &one, NULL, &status);
		fits_update_key(fptr, TDOUBLE, "BZERO", &zero, NULL, &status);
		fits_update_key(fptr, TLOGICAL, "LSF", &t, "Living Stack Format file", &status);
		fits_update_key(fptr, TSTRING, "LSFVER", (char *) LSF_VERSION, "LSF FITS profile version", &status);
		write_identity(fptr, NULL, uid, gen, &status);
		fits_update_key(fptr, TSTRING, "LSFSENS", (char *) (s->is_cfa ? "CFA" : "MONO"), "LSF sensor class", &status);
		if (s->is_cfa)
			fits_update_key(fptr, TSTRING, "LSFBAYR", (void *) s->bayer, "LSF Bayer pattern of the sensor", &status);
		else delete_key_quiet(fptr, "LSFBAYR");
		fits_update_key(fptr, TSTRING, "LSFDMOD", (char *) (s->is_cfa ? "DEBAYER" : "NONE"), "LSF demosaic mode", &status);
		fits_update_key(fptr, TLOGICAL, "LSFNORM", s->outnorm ? &t : &f,
				"Stack is output-normalized", &status);
		if (s->outnorm) {
			/* exact round trip: 17 significant digits */
			fits_update_key_dbl(fptr, "LSFNMIN", s->outnorm_min, -17, "output normalization min", &status);
			fits_update_key_dbl(fptr, "LSFNMAX", s->outnorm_max, -17, "output normalization max", &status);
		} else {
			delete_key_quiet(fptr, "LSFNMIN");
			delete_key_quiet(fptr, "LSFNMAX");
		}
		delete_key_quiet(fptr, "LSFBUSY");
		gchar *hash = lsf_hash_float_data(primary->fdata, npix * nchan);
		fits_update_key(fptr, TSTRING, "LSFHASH", hash, "SHA-256 of the primary data", &status);
		g_free(hash);
	}
	if (status) {
		report_fits_error(status);
		goto end;
	}
	wcs = collect_wcs_cards(fptr);

	/* LSF_ACC */
	long anaxes[3] = { h->rx, h->ry, LSF_NB_PL * nchan };
	fits_create_img(fptr, FLOAT_IMG, 3, anaxes, &status);
	write_identity(fptr, LSF_EXT_ACC, uid, gen, &status);
	for (int k = 0, p = 1; k < LSF_NB_PL; k++)
		for (int c = 0; c < nchan; c++, p++) {
			char key[FLEN_KEYWORD];
			gchar *name = plane_name(k, c, nchan);
			g_snprintf(key, sizeof key, "LSFPL%d", p);
			fits_update_key(fptr, TSTRING, key, name, NULL, &status);
			g_free(name);
		}
	write_cards(fptr, wcs, &status);
	for (int k = 0, p = 1; k < LSF_NB_PL && !status; k++)
		for (int c = 0; c < nchan && !status; c++, p++) {
			long fpixel[3] = { 1, 1, p };
			fits_write_pix(fptr, TFLOAT, fpixel, npix, h->plane[k][c], &status);
		}

	/* LSF_COUNT */
	long cnaxes[3] = { h->rx, h->ry, nchan };
	fits_create_img(fptr, USHORT_IMG, naxis, cnaxes, &status);
	write_identity(fptr, LSF_EXT_COUNT, uid, gen, &status);
	if (nchan == 3) {
		fits_update_key(fptr, TSTRING, "LSFPL1", (char *) "COUNT_R", NULL, &status);
		fits_update_key(fptr, TSTRING, "LSFPL2", (char *) "COUNT_G", NULL, &status);
		fits_update_key(fptr, TSTRING, "LSFPL3", (char *) "COUNT_B", NULL, &status);
	}
	if (h->count_saturated) {
		int t = TRUE;
		fits_update_key(fptr, TLOGICAL, "LSFCSAT", &t, "a count saturated at 65535", &status);
	}
	write_cards(fptr, wcs, &status);
	for (int c = 0; c < nchan && !status; c++) {
		long fpixel[3] = { 1, 1, c + 1 };
		fits_write_pix(fptr, TUSHORT, fpixel, npix, h->count[c], &status);
	}

	/* LSF_SETTINGS */
	write_settings_hdu(fptr, s, nchan, uid, gen, &status);
	if (status) {
		report_fits_error(status);
		goto end;
	}

	/* unknown HDUs are preserved, before the ledger (C7, Section 3) */
	if (old) {
		if (fits_open_diskfile(&oldf, old->path, READONLY, &status)) {
			report_fits_error(status);
			goto end;
		}
		for (guint i = 0; i < old->unknown_hdus->len && !status; i++) {
			fits_movabs_hdu(oldf, g_array_index(old->unknown_hdus, int, i), NULL, &status);
			fits_copy_hdu(oldf, fptr, 0, &status);
		}
	}

	/* LSF_LEDGER, last */
	long first_row = 1;
	if (old) {
		fits_movabs_hdu(oldf, old->ledger_hdu, NULL, &status);
		fits_copy_hdu(oldf, fptr, 0, &status);	// keeps unknown columns (C7)
		first_row = old->ledger_rows + 1;
	} else {
		create_ledger_hdu(fptr, nchan, &status);
	}
	write_identity(fptr, NULL, uid, gen, &status);
	write_ledger_rows(fptr, new_rows, first_row, nchan, batch, gen, added, &status);
	if (status) {
		report_fits_error(status);
		goto end;
	}

	/* DATASUM and CHECKSUM on every HDU */
	int hdunum = 0;
	fits_get_num_hdus(fptr, &hdunum, &status);
	for (int i = 1; i <= hdunum && !status; i++) {
		fits_movabs_hdu(fptr, i, NULL, &status);
		fits_write_chksum(fptr, &status);
	}
	if (status) {
		report_fits_error(status);
		goto end;
	}
	fits_close_file(fptr, &status);
	fptr = NULL;
	primary->fptr = NULL;
	if (oldf) {
		int st = 0;
		fits_close_file(oldf, &st);
		oldf = NULL;
	}
	if (status) {
		report_fits_error(status);
		goto end;
	}
	if (flush_to_disk(tmp))
		siril_log_warning(_("LSF: could not flush %s to disk\n"), tmp);

	/* atomic replace, keeping the temporary file on failure */
	if (backup && g_file_test(path, G_FILE_TEST_EXISTS)) {
		gchar *bak = g_strdup_printf("%s.bak", path);
		g_unlink(bak);
		if (g_rename(path, bak))
			siril_log_warning(_("LSF: could not keep a backup copy as %s\n"), bak);
		else siril_log_message(_("LSF: previous file kept as %s\n"), bak);
		g_free(bak);
	}
	if (replace_file(tmp, path)) {
		siril_log_error(_("LSF: could not replace %s (is it open in another program?). "
					"The new file was kept as %s\n"), path, tmp);
		g_free(tmp);
		tmp = NULL;
		goto end;
	}
	retval = 0;

end:
	if (fptr) {
		int st = 0;
		fits_close_file(fptr, &st);
		primary->fptr = NULL;
		if (tmp) g_unlink(tmp);
	}
	if (oldf) {
		int st = 0;
		fits_close_file(oldf, &st);
	}
	if (wcs)
		g_ptr_array_free(wcs, TRUE);
	g_free(tmp);
	g_free(added);
	return retval;
}

/* removes the LSF keywords that Siril kept as unknown keywords of a header
 * copied from an LSF primary, so a derived image does not claim to be LSF */
void lsf_strip_lsf_keys(fits *fit) {
	if (!fit->unknown_keys)
		return;
	gchar **cards = g_strsplit(fit->unknown_keys, "\n", -1);
	GString *kept = g_string_new(NULL);
	for (int i = 0; cards[i]; i++) {
		char key[FLEN_KEYWORD];
		int len = 0, status = 0;
		if (!cards[i][0])
			continue;
		fits_get_keyname(cards[i], key, &len, &status);
		gboolean lsf = FALSE;
		for (int k = 0; !status && k < (int) G_N_ELEMENTS(lsf_keys); k++)
			lsf |= !strcmp(key, lsf_keys[k]);
		if (lsf)
			continue;
		g_string_append(kept, cards[i]);
		g_string_append_c(kept, '\n');
	}
	g_strfreev(cards);
	g_free(fit->unknown_keys);
	fit->unknown_keys = g_string_free(kept, FALSE);
}

/* Flattened master (Section 10.9): HDU 0 without the LSF keywords */
int lsf_export_master(const char *lsf_path, const char *out_path) {
	fitsfile *in = NULL, *out = NULL;
	int status = 0, gen = 0;
	char uid[FLEN_VALUE] = { 0 };
	if (fits_open_diskfile(&in, lsf_path, READONLY, &status)) {
		report_fits_error(status);
		return 1;
	}
	read_uid_gen(in, uid, &gen);
	g_unlink(out_path);
	if (siril_fits_create_diskfile(&out, out_path, &status)) {
		report_fits_error(status);
		fits_close_file(in, &status);
		return 1;
	}
	fits_copy_hdu(in, out, 0, &status);
	for (int i = 0; i < (int) G_N_ELEMENTS(lsf_keys); i++)
		delete_key_quiet(out, lsf_keys[i]);
	gchar *hist = g_strdup_printf("Exported from LSF file %s, generation %d", uid, gen);
	fits_write_history(out, hist, &status);
	g_free(hist);
	fits_write_chksum(out, &status);
	int st = 0;
	fits_close_file(out, &status);
	fits_close_file(in, &st);
	if (status) {
		report_fits_error(status);
		return 1;
	}
	return 0;
}
