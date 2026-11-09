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

/* LSF_SETTINGS rows (FITS profile v0.4, Section 6): KEY, TYPE, VALUE */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>

#include "core/siril_log.h"
#include "lsf/lsf.h"
#include "lsf/lsf_settings.h"

#include "yyjson.h"

static const char *rej_names[] = { "none", "percentile", "sigma", "mad", "sigmedian",
	"winsorized", "linearfit", "gesdt" };
static const char *norm_names[] = { "none", "add", "mul", "addscale", "mulscale" };
static const char *weight_names[] = { "none", "nbstars", "wfwhm", "noise", "nbstack" };
static const char *transf_names[] = { "shift", "similarity", "affine", "homography" };
static const char *interp_names[] = { "nearest", "linear", "cubic", "area", "lanczos4", "none" };

const char *lsf_rejection_name(rejection r) { return rej_names[r]; }
gboolean lsf_rejection_from_name(const char *name, rejection *r) {
	for (int i = 0; i < (int) G_N_ELEMENTS(rej_names); i++)
		if (!g_ascii_strcasecmp(rej_names[i], name)) { *r = (rejection) i; return TRUE; }
	return FALSE;
}
const char *lsf_norm_name(normalization n) { return norm_names[n]; }
gboolean lsf_norm_from_name(const char *name, normalization *n) {
	for (int i = 0; i < (int) G_N_ELEMENTS(norm_names); i++)
		if (!g_ascii_strcasecmp(norm_names[i], name)) { *n = (normalization) i; return TRUE; }
	return FALSE;
}
/* align.disto: the distortion model of the appends, a file of Siril's
 * register -disto=file, the master pattern of -disto=master, or the frames'
 * own SIP terms of an astrometric registration (JSON: source and file, empty
 * for the frames) */
gchar *lsf_disto_to_json(const disto_params *d) {
	yyjson_mut_doc *doc = yyjson_mut_doc_new(NULL);
	yyjson_mut_val *root = yyjson_mut_obj(doc);
	yyjson_mut_doc_set_root(doc, root);
	yyjson_mut_obj_add_str(doc, root, "source", d->index == DISTO_MASTER ? "master" : d->index == DISTO_FILES ? "frames" : "file");
	yyjson_mut_obj_add_str(doc, root, "file", d->filename ? d->filename : "");
	char *json = yyjson_mut_write(doc, 0, NULL);
	gchar *out = g_strdup(json ? json : "{}");
	free(json);
	yyjson_mut_doc_free(doc);
	return out;
}

int lsf_disto_from_json(const char *json, disto_params *d) {
	yyjson_doc *doc = yyjson_read(json, strlen(json), 0);
	if (!doc) return 1;
	yyjson_val *root = yyjson_doc_get_root(doc);
	const char *src = yyjson_get_str(yyjson_obj_get(root, "source"));
	const char *file = yyjson_get_str(yyjson_obj_get(root, "file"));
	int ret = !src || !file || (strcmp(src, "master") && strcmp(src, "file") && strcmp(src, "frames"));
	if (!ret) {
		g_free(d->filename);
		d->index = !strcmp(src, "master") ? DISTO_MASTER : !strcmp(src, "frames") ? DISTO_FILES : DISTO_FILE;
		d->filename = d->index == DISTO_FILES ? NULL : g_strdup(file);
	}
	yyjson_doc_free(doc);
	return ret;
}

void lsf_disto_copy(disto_params *dst, const disto_params *src) {
	g_free(dst->filename);
	dst->index = src->index;
	dst->filename = g_strdup(src->filename);
	dst->velocity = src->velocity;
}

/* how a model is shown, and given back to lsf_append (-disto=) */
gchar *lsf_disto_describe(const disto_params *d) {
	if (d->index == DISTO_MASTER)
		return g_strdup_printf("master (%s)", d->filename ? d->filename : "");
	if (d->index == DISTO_FILE)
		return g_strdup_printf("file %s", d->filename ? d->filename : "");
	if (d->index == DISTO_FILES)
		return g_strdup(_("the SIP terms of each frame's own WCS"));
	return g_strdup("none");
}

gboolean lsf_disto_same(const disto_params *a, const disto_params *b) {
	return a->index == b->index && !g_strcmp0(a->filename, b->filename);
}

static const char *weight_model_names[] = { "range", "inverse" };
const char *lsf_weight_model_name(int m) { return weight_model_names[m]; }
gboolean lsf_weight_model_from_name(const char *name, int *m) {
	for (int i = 0; i < (int) G_N_ELEMENTS(weight_model_names); i++)
		if (!g_ascii_strcasecmp(weight_model_names[i], name)) { *m = i; return TRUE; }
	return FALSE;
}

const char *lsf_weight_name(weightingType w) { return weight_names[w]; }
gboolean lsf_weight_from_name(const char *name, weightingType *w) {
	for (int i = 0; i < (int) G_N_ELEMENTS(weight_names); i++)
		if (!g_ascii_strcasecmp(weight_names[i], name)) { *w = (weightingType) i; return TRUE; }
	return FALSE;
}

/* methods whose thresholds are expressed in sigma units, Section 10.6 */
gboolean lsf_rejection_has_sigma(rejection r) {
	return r == SIGMA || r == MAD || r == SIGMEDIAN || r == WINSORIZED || r == LINEARFIT;
}

void lsf_setting_row_free(gpointer p) {
	lsf_setting_row *r = p;
	if (!r) return;
	g_free(r->key);
	g_free(r->type);
	g_free(r->value);
	g_free(r);
}

void lsf_settings_init(lsf_settings *s) {
	memset(s, 0, sizeof(lsf_settings));
	g_strlcpy(s->bayer, "NONE", sizeof s->bayer);
	s->binning = 1;
	s->oversample = 1;
	s->transform = HOMOGRAPHY_TRANSFORMATION;
	s->interpolation = OPENCV_LANCZOS4;
	s->min_pairs = 10;
	s->rej = WINSORIZED;
	s->sig[0] = s->sig[1] = 3.0;
	s->hist_low = s->hist_high = 3.0;
	s->hist_margin = 1.0;
	s->weighting = NO_WEIGHT;
	s->fwhm_norm = s->stars_norm = 1.0;
	s->weight_model = LSF_WMODEL_RANGE;
	s->weight_floor = LSF_WEIGHT_FLOOR;
	s->weight_range = 1.0;
	s->norm = ADDITIVE_SCALING;
	s->unknown = g_ptr_array_new_with_free_func(lsf_setting_row_free);
}

void lsf_settings_free(lsf_settings *s) {
	g_free(s->findstar_json);
	g_free(s->filter);
	g_free(s->created_utc);
	g_free(s->disto.filename);
	if (s->unknown)
		g_ptr_array_free(s->unknown, TRUE);
	s->findstar_json = s->filter = s->created_utc = s->disto.filename = NULL;
	s->unknown = NULL;
}

static void add_row(GPtrArray *rows, const char *key, const char *type, gchar *value) {
	lsf_setting_row *r = g_new(lsf_setting_row, 1);
	r->key = g_strdup(key);
	r->type = g_strdup(type);
	r->value = value;	// owned
	g_ptr_array_add(rows, r);
}

static gchar *fmt_double(double v) {
	gchar buf[G_ASCII_DTOSTR_BUF_SIZE];
	return g_strdup(g_ascii_formatd(buf, sizeof buf, "%.17g", v));
}

static void add_str(GPtrArray *rows, const char *key, const char *v) { add_row(rows, key, "STR", g_strdup(v ? v : "")); }
static void add_int(GPtrArray *rows, const char *key, int v) { add_row(rows, key, "INT", g_strdup_printf("%d", v)); }
static void add_flt(GPtrArray *rows, const char *key, double v) { add_row(rows, key, "FLT", fmt_double(v)); }
static void add_bool(GPtrArray *rows, const char *key, gboolean v) { add_row(rows, key, "BOOL", g_strdup(v ? "T" : "F")); }

static gchar *json_array(const double *v, int n) {
	GString *s = g_string_new("[");
	for (int i = 0; i < n; i++) {
		gchar *d = fmt_double(v[i]);
		g_string_append_printf(s, "%s%s", i ? "," : "", d);
		g_free(d);
	}
	g_string_append_c(s, ']');
	return g_string_free(s, FALSE);
}

GPtrArray *lsf_settings_to_rows(const lsf_settings *s, int nchan) {
	GPtrArray *rows = g_ptr_array_new_with_free_func(lsf_setting_row_free);
	add_str(rows, "sensor.class", s->is_cfa ? "CFA" : "MONO");
	add_str(rows, "sensor.bayer", s->bayer);
	add_int(rows, "sensor.binning", s->binning);
	add_str(rows, "cfa.mode", !s->is_cfa ? "NONE" : s->drizzled ? "CFADRIZ" : "DEBAYER");
	/* written back as read, so a rewrite never changes what the grid is.
	 * The drizzle pixel fraction and kernel are choices of each append, not
	 * stored (Section 6) */
	add_int(rows, "grid.oversample", s->oversample);
	add_int(rows, "grid.origin_x", s->origin_x);
	add_int(rows, "grid.origin_y", s->origin_y);
	add_bool(rows, "grid.undistorted", s->undistorted);
	add_str(rows, "align.transform", transf_names[s->transform < 0 ? HOMOGRAPHY_TRANSFORMATION : s->transform]);
	add_str(rows, "align.interpolation", interp_names[s->interpolation]);
	add_int(rows, "align.min_pairs", s->min_pairs);
	add_row(rows, "align.findstar", "JSON", g_strdup(s->findstar_json ? s->findstar_json : "{}"));
	if (s->undistorted && s->disto.index != DISTO_UNDEF)
		add_row(rows, "align.disto", "JSON", lsf_disto_to_json(&s->disto));
	add_str(rows, "reject.method", lsf_rejection_name(s->rej));
	gchar *lo = fmt_double(s->sig[0]), *hi = fmt_double(s->sig[1]);
	add_row(rows, "reject.params", "JSON", g_strdup_printf("{\"low\":%s,\"high\":%s}", lo, hi));
	g_free(lo); g_free(hi);
	add_flt(rows, "reject.hist_low", s->hist_low);
	add_flt(rows, "reject.hist_high", s->hist_high);
	add_flt(rows, "reject.hist_margin", s->hist_margin);
	add_str(rows, "weight.scheme", lsf_weight_name(s->weighting));
	switch (s->weighting) {
	case NOISE_WEIGHT:
		add_row(rows, "weight.noise_norm", "JSON", json_array(s->noise_norm, nchan));
		break;
	case WFWHM_WEIGHT:
		add_flt(rows, "weight.fwhm_min", s->fwhm_min);
		add_flt(rows, "weight.fwhm_max", s->fwhm_max);
		add_flt(rows, "weight.fwhm_norm", s->fwhm_norm);
		break;
	case NBSTARS_WEIGHT:
		add_flt(rows, "weight.stars_min", s->stars_min);
		add_flt(rows, "weight.stars_max", s->stars_max);
		add_flt(rows, "weight.stars_norm", s->stars_norm);
		break;
	default:
		break;
	}
	/* the model, floor and span multiplier of the schemes that freeze a range */
	if (s->weighting == WFWHM_WEIGHT || s->weighting == NBSTARS_WEIGHT) {
		add_str(rows, "weight.model", lsf_weight_model_name(s->weight_model));
		add_flt(rows, "weight.floor", s->weight_floor);
		add_flt(rows, "weight.range", s->weight_range);
	}
	if (s->drizzled)
		add_row(rows, "weight.driz_norm", "JSON", json_array(s->driz_norm, nchan));
	if (s->drizzled && s->driz_flat)	// absent means no flat, as in files made before it
		add_bool(rows, "weight.driz_flat", TRUE);
	add_str(rows, "norm.method", lsf_norm_name(s->norm));
	add_str(rows, "norm.estimator", s->lite_norm ? "lite" : "ikss");
	add_bool(rows, "norm.rgb_equal", s->rgb_equal);
	if (s->norm != NO_NORM) {
		gchar *l = json_array(s->base_loc, nchan), *sc = json_array(s->base_scale, nchan);
		add_row(rows, "norm.baseline", "JSON", g_strdup_printf("{\"loc\":%s,\"scale\":%s}", l, sc));
		g_free(l); g_free(sc);
	}
	add_bool(rows, "outnorm.enabled", s->outnorm);
	if (s->outnorm) {
		add_flt(rows, "outnorm.min", s->outnorm_min);
		add_flt(rows, "outnorm.max", s->outnorm_max);
	}
	add_bool(rows, "rejmaps.enabled", s->rejmaps);
	add_str(rows, "filter", s->filter);
	add_str(rows, "created_utc", s->created_utc);
	add_str(rows, "writer", "Siril " PACKAGE_VERSION);
	/* C7: unknown rows, unless this write set the same key */
	if (s->unknown) {
		for (guint i = 0; i < s->unknown->len; i++) {
			lsf_setting_row *u = g_ptr_array_index(s->unknown, i);
			gboolean set = FALSE;
			for (guint j = 0; j < rows->len && !set; j++)
				set = !strcmp(((lsf_setting_row *) g_ptr_array_index(rows, j))->key, u->key);
			if (!set)
				add_row(rows, u->key, u->type, g_strdup(u->value));
		}
	}
	return rows;
}

static gboolean parse_bool(const char *v) {
	return v && (v[0] == 'T' || v[0] == 't' || v[0] == '1');
}

static int parse_json_array(const char *json, double *out, int n) {
	yyjson_doc *doc = yyjson_read(json, strlen(json), 0);
	if (!doc) return 1;
	yyjson_val *arr = yyjson_doc_get_root(doc);
	int ret = (!yyjson_is_arr(arr) || (int) yyjson_arr_size(arr) < n);
	for (int i = 0; !ret && i < n; i++)
		out[i] = yyjson_get_num(yyjson_arr_get(arr, i));
	yyjson_doc_free(doc);
	return ret;
}

static int parse_baseline(const char *json, lsf_settings *s, int n) {
	yyjson_doc *doc = yyjson_read(json, strlen(json), 0);
	if (!doc) return 1;
	yyjson_val *root = yyjson_doc_get_root(doc);
	yyjson_val *loc = yyjson_obj_get(root, "loc"), *sc = yyjson_obj_get(root, "scale");
	int ret = (!yyjson_is_arr(loc) || !yyjson_is_arr(sc) ||
			(int) yyjson_arr_size(loc) < n || (int) yyjson_arr_size(sc) < n);
	for (int i = 0; !ret && i < n; i++) {
		s->base_loc[i] = yyjson_get_num(yyjson_arr_get(loc, i));
		s->base_scale[i] = yyjson_get_num(yyjson_arr_get(sc, i));
	}
	yyjson_doc_free(doc);
	return ret;
}

static int parse_reject_params(const char *json, lsf_settings *s) {
	yyjson_doc *doc = yyjson_read(json, strlen(json), 0);
	if (!doc) return 1;
	yyjson_val *root = yyjson_doc_get_root(doc);
	yyjson_val *lo = yyjson_obj_get(root, "low"), *hi = yyjson_obj_get(root, "high");
	int ret = !yyjson_is_num(lo) || !yyjson_is_num(hi);
	if (!ret) {
		s->sig[0] = yyjson_get_num(lo);
		s->sig[1] = yyjson_get_num(hi);
	}
	yyjson_doc_free(doc);
	return ret;
}

static int name_index(const char **names, int n, const char *v) {
	for (int i = 0; i < n; i++)
		if (!g_ascii_strcasecmp(names[i], v)) return i;
	return -1;
}

/* fills s (already initialized) from the rows read in the file; unknown
 * keys are kept. Returns non-zero on a malformed required value. */
int lsf_settings_from_rows(lsf_settings *s, GPtrArray *rows, int nchan) {
	int err = 0;
	for (guint i = 0; i < rows->len; i++) {
		lsf_setting_row *r = g_ptr_array_index(rows, i);
		const char *k = r->key, *v = r->value ? r->value : "";
		int idx;
		if (!strcmp(k, "sensor.class")) s->is_cfa = !g_ascii_strcasecmp(v, "CFA");
		else if (!strcmp(k, "sensor.bayer")) g_strlcpy(s->bayer, v, sizeof s->bayer);
		else if (!strcmp(k, "sensor.binning")) s->binning = atoi(v);
		else if (!strcmp(k, "grid.origin_x")) s->origin_x = atoi(v);
		else if (!strcmp(k, "grid.origin_y")) s->origin_y = atoi(v);
		else if (!strcmp(k, "grid.oversample")) {
			gchar *end;
			gint64 o = g_ascii_strtoll(v, &end, 10);
			if (end == v || o < 1 || o > 16) err = 1; else s->oversample = (int) o;
		}
		else if (!strcmp(k, "grid.expand") || !strcmp(k, "grid.expand_step") || !strcmp(k, "grid.pixfrac")
				|| !strcmp(k, "grid.kernel") || !strcmp(k, "grid.reference")) {
			/* written by earlier versions, dropped (Section 6): chosen per append */
		}
		else if (!strcmp(k, "align.transform")) {
			idx = name_index(transf_names, G_N_ELEMENTS(transf_names), v);
			if (idx < 0) err = 1; else s->transform = (transformation_type) idx;
		}
		else if (!strcmp(k, "align.interpolation")) {
			idx = name_index(interp_names, G_N_ELEMENTS(interp_names), v);
			if (idx < 0) err = 1; else s->interpolation = (opencv_interpolation) idx;
		}
		else if (!strcmp(k, "align.min_pairs")) s->min_pairs = atoi(v);
		else if (!strcmp(k, "align.findstar")) { g_free(s->findstar_json); s->findstar_json = g_strdup(v); }
		else if (!strcmp(k, "grid.undistorted")) s->undistorted = parse_bool(v);
		else if (!strcmp(k, "align.disto")) err |= lsf_disto_from_json(v, &s->disto);
		else if (!strcmp(k, "reject.method")) err |= !lsf_rejection_from_name(v, &s->rej);
		else if (!strcmp(k, "reject.params")) err |= parse_reject_params(v, s);
		else if (!strcmp(k, "reject.hist_low")) s->hist_low = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "reject.hist_high")) s->hist_high = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "reject.hist_margin")) s->hist_margin = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "weight.scheme")) err |= !lsf_weight_from_name(v, &s->weighting);
		else if (!strcmp(k, "weight.noise_norm")) err |= parse_json_array(v, s->noise_norm, nchan);
		else if (!strcmp(k, "weight.driz_norm")) {	// a drizzled file, CFA or mono
			err |= parse_json_array(v, s->driz_norm, nchan);
			s->drizzled = TRUE;
		}
		else if (!strcmp(k, "weight.driz_flat")) s->driz_flat = parse_bool(v);
		else if (!strcmp(k, "weight.fwhm_min")) s->fwhm_min = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "weight.fwhm_max")) s->fwhm_max = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "weight.fwhm_norm")) s->fwhm_norm = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "weight.model")) err |= !lsf_weight_model_from_name(v, &s->weight_model);
		else if (!strcmp(k, "weight.floor")) s->weight_floor = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "weight.range")) s->weight_range = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "weight.stars_min")) s->stars_min = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "weight.stars_max")) s->stars_max = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "weight.stars_norm")) s->stars_norm = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "norm.method")) err |= !lsf_norm_from_name(v, &s->norm);
		else if (!strcmp(k, "norm.estimator")) s->lite_norm = !g_ascii_strcasecmp(v, "lite");
		else if (!strcmp(k, "norm.rgb_equal")) s->rgb_equal = parse_bool(v);
		else if (!strcmp(k, "norm.baseline")) err |= parse_baseline(v, s, nchan);
		else if (!strcmp(k, "outnorm.enabled")) s->outnorm = parse_bool(v);
		else if (!strcmp(k, "outnorm.min")) s->outnorm_min = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "outnorm.max")) s->outnorm_max = g_ascii_strtod(v, NULL);
		else if (!strcmp(k, "rejmaps.enabled")) s->rejmaps = parse_bool(v);
		else if (!strcmp(k, "filter")) { g_free(s->filter); s->filter = g_strdup(v); }
		else if (!strcmp(k, "created_utc")) { g_free(s->created_utc); s->created_utc = g_strdup(v); }
		else if (!strcmp(k, "cfa.mode")) s->drizzled |= !g_ascii_strcasecmp(v, "CFADRIZ");
		else if (!strcmp(k, "writer")) {
			/* a fixed value in this version, rewritten on save */
		}
		else {
			lsf_setting_row *u = g_new(lsf_setting_row, 1);
			u->key = g_strdup(k);
			u->type = g_strdup(r->type);
			u->value = g_strdup(v);
			g_ptr_array_add(s->unknown, u);
		}
	}
	/* the history test thresholds are in sigma units: 0 would reject every
	 * sample on that side of the mean */
	if (!(s->hist_low > 0.0) || !(s->hist_high > 0.0) || !(s->hist_margin >= 1.0)) {
		siril_log_error(_("LSF: the history test settings are not valid (reject.hist_low and reject.hist_high "
					"must be above 0, reject.hist_margin at least 1)\n"));
		err = 1;
	}
	if (!(s->weight_floor > 0.0 && s->weight_floor < 1.0) || !(s->weight_range >= 1.0)) {
		siril_log_error(_("LSF: the weight settings are not valid (weight.floor must be between 0 and 1, "
					"weight.range at least 1)\n"));
		err = 1;
	}
	for (int c = 0; s->drizzled && c < nchan; c++)
		if (!(s->driz_norm[c] > 0.0)) {
			siril_log_error(_("LSF: a drizzled file needs weight.driz_norm, above 0 for every channel\n"));
			err = 1;
			break;
		}
	return err;
}

/* applies the star detection parameters stored in align.findstar to sf */
void lsf_findstar_from_json(const char *json, star_finder_params *sf) {
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
