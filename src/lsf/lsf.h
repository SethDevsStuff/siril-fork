#ifndef SRC_LSF_LSF_H_
#define SRC_LSF_LSF_H_

/* Living Stack Format (LSF), FITS profile v0.4.
 * See lsf_poc/live stack fits.md for the specification. */

#include <glib.h>
#include "core/siril.h"
#include "core/settings.h"
#include "stacking/stacking.h"
#include "registration/registration.h"

#define LSF_VERSION "0.4"
#define LSF_EXT_ACC "LSF_ACC"
#define LSF_EXT_COUNT "LSF_COUNT"
#define LSF_EXT_REJ "LSF_REJ"
#define LSF_EXT_SETTINGS "LSF_SETTINGS"
#define LSF_EXT_LEDGER "LSF_LEDGER"
#define LSF_COUNT_MAX 65535
#define LSF_LOW_STACK_WARN 25
#define LSF_NEFF_MIN (1.0 + 1e-4)
/* least wfwhm or nbstars weight before normalization (Section 10.3), 1 being
 * the best frame of the first batch: a frame at or beyond its worst still counts */
#define LSF_WEIGHT_FLOOR 0.01

/* tile compression of the extension images (Section 5.3): LSF_ACC with
 * GZIP_1 or GZIP_2 and no quantization, the count planes with RICE_1 */
typedef enum { LSF_COMP_NONE, LSF_COMP_GZIP1, LSF_COMP_GZIP2 } lsf_compression;

/* planes of the accumulator, in LSF_ACC order */
enum { LSF_PL_MEAN, LSF_PL_W, LSF_PL_M2, LSF_PL_W2, LSF_NB_PL };

/* accumulator history: one plane per (kind, channel), FITS memory order
 * (the same order as fits.fpdata) */
typedef struct {
	int rx, ry, nchan;
	float *plane[LSF_NB_PL][3];
	WORD *count[3];
	gboolean count_saturated;
	/* rejection maps (Section 5.4), NULL when the file has no LSF_REJ:
	 * rejected samples per pixel, low then high */
	WORD *rej[2][3];
	gboolean rej_saturated;
	int rej_gen;			// LSFRGEN, the first generation counted
} lsf_history;

/* one row of LSF_SETTINGS that this writer does not know, kept (C7) */
typedef struct {
	gchar *key, *type, *value;
} lsf_setting_row;

typedef struct {
	/* sensor */
	gboolean is_cfa;		// sensor.class
	gchar bayer[16];		// sensor.bayer, "NONE" for mono
	int binning;			// sensor.binning
	/* grid */
	int origin_x, origin_y;		// grid.origin_*, total padding before the data (Section 10.5)
	int oversample;			// grid.oversample, the drizzle scale: this writer merges onto 1 only
	/* registration */
	transformation_type transform;	// align.transform
	opencv_interpolation interpolation;	// align.interpolation
	int min_pairs;			// align.min_pairs
	gchar *findstar_json;		// align.findstar
	/* rejection */
	rejection rej;			// reject.method
	double sig[2];			// reject.params
	double hist_low, hist_high, hist_margin;
	/* weighting, with constants frozen from the first batch (§10.3) */
	weightingType weighting;	// weight.scheme
	double noise_norm[3];		// per channel
	double fwhm_min, fwhm_max, fwhm_norm;
	double stars_min, stars_max, stars_norm;
	/* normalization against a frozen baseline frame */
	normalization norm;		// norm.method
	gboolean lite_norm;		// norm.estimator: lite or ikss
	gboolean rgb_equal;		// norm.rgb_equal
	double base_loc[3], base_scale[3];	// norm.baseline
	/* output normalization */
	gboolean outnorm;
	double outnorm_min, outnorm_max;
	gboolean rejmaps;
	gchar *filter;
	gchar *created_utc;
	GPtrArray *unknown;		// of lsf_setting_row*
} lsf_settings;

/* values from the existing file needed to rewrite it */
typedef struct {
	gchar *path;
	gchar uid[64];
	int gen;
	int nchan, rx, ry;
	lsf_settings settings;
	lsf_history hist;
	fits primary;			// metadata only (no pixels), for the header
	int ledger_hdu;			// 1-based HDU number of LSF_LEDGER
	GArray *unknown_hdus;		// of int, 1-based HDU numbers to preserve (C7)
	long ledger_rows;
	int max_ledger_gen;
	int nb_batches;
	double outnorm_min, outnorm_max;	// LSFNMIN/LSFNMAX
	gboolean stack_modified;	// check 8 failed (only when tolerated)
	float *stack_data;		// HDU 0 pixels when stack_modified, file order
	double xfrm_shift[2];		// grid expansion of this commit: shift of the REG_XFRM already in the ledger
	lsf_compression compression;	// of LSF_ACC, as found in the file
	gchar acc_codec[FLEN_VALUE], count_codec[FLEN_VALUE], rej_codec[FLEN_VALUE];	// ZCMPTYPE, "none" when not compressed
	long long acc_bytes, count_bytes, rej_bytes;	// size of the data units as stored
} lsf_file;

/* one frame of a batch, for the ledger */
typedef struct {
	gchar *frame_id;
	gchar *filename;
	gchar *date_obs;	// ISO 8601, may be NULL
	double exptime;
	int gain, offset;	// INT_MIN when unknown
	int xbin, ybin;
	float ccd_temp;
	gchar *filter;
	gchar *orighdr;
	gchar reg_type[9];
	double xfrm[9];
	gchar *align_json;
	double norm_sc[3], norm_off[3];
	gchar norm_mth[3][13];
	double weight[3];
	float fwhm, wfwhm, round, regq, regbg;
	int nstars;
	float bkg_med[3], noise[3], noise_n[3];
	float st_mean[3], st_sig[3], st_mad[3], st_sqbw[3], st_loc[3], st_scl[3];
} lsf_ledger_row;

/* sides of a grid expansion, in FITS storage order (Section 10.5): FIRST
 * pads before FITS row 1, LAST after row NAXIS2. Registration coordinates
 * run the other way in y, so a frame above y = 0 there needs LAST rows. */
enum { LSF_PAD_LEFT, LSF_PAD_RIGHT, LSF_PAD_FIRST, LSF_PAD_LAST };

/* lsf_math.c: pure functions, unit-tested */
void lsf_batch_stats(const float *x, const double *w, int n,
		double *mean, double *W, double *M2, double *W2);
void lsf_merge(float *mean_a, float *w_a, float *m2_a, float *w2_a, WORD *c_a,
		double mean_b, double w_b, double m2_b, double w2_b, int c_b, gboolean *saturated);
void lsf_count_add(WORD *c, int n, gboolean *saturated);
double lsf_clip_correction(double kl, double kh);
double lsf_sigma_pred2(weightingType scheme, double W, double M2, double W2, int N, double w_new);
double lsf_sigma_hist(double sigma_pred2, double c, double sigma_frame, double margin);
gboolean lsf_outnorm_range(float *const mean[3], float *const W[3], size_t n, int nchan,
		double *mini, double *maxi);
void lsf_pixel_stats(weightingType scheme, double W, double M2, double W2, int N,
		double *neff, double *sigma2, double *var_mean);
float lsf_stack_value(float mean, float w, gboolean norm, double mini, double range);
size_t lsf_adopt_plane(float *mean, const float *w, const float *disp, size_t n,
		gboolean norm, double mini, double maxi, size_t *nonfinite);
void lsf_footprint_init(double box[4]);
void lsf_footprint_add(const double h[9], int fx, int fy, double box[4]);
gboolean lsf_grid_overhang(const double box[4], int rx, int ry, int need[4]);
void lsf_expand_pads(const int need[4], int step, int pad[4]);
void lsf_xfrm_shift(double h[9], double dx, double dy);
void lsf_pad_plane(const void *src, void *dst, size_t elsize, int rx, int ry, const int pad[4]);

/* lsf_weights.c */
int lsf_freeze_weights(lsf_settings *s, int nb_frames, int nchan, const double *noise_raw,
		const float *wfwhm, const int *nstars);
double lsf_frame_weight(const lsf_settings *s, int chan, double noise_raw, float wfwhm, int nstars);
double lsf_noise_raw_weight(double pscale, double bgnoise);

/* lsf_settings.c */
void lsf_settings_init(lsf_settings *s);
void lsf_settings_free(lsf_settings *s);
const char *lsf_rejection_name(rejection r);
gboolean lsf_rejection_from_name(const char *name, rejection *r);
const char *lsf_norm_name(normalization n);
gboolean lsf_norm_from_name(const char *name, normalization *n);
const char *lsf_weight_name(weightingType w);
gboolean lsf_weight_from_name(const char *name, weightingType *w);
gboolean lsf_rejection_has_sigma(rejection r);
GPtrArray *lsf_settings_to_rows(const lsf_settings *s, int nchan);
int lsf_settings_from_rows(lsf_settings *s, GPtrArray *rows, int nchan);
void lsf_setting_row_free(gpointer p);

/* lsf_fits.c */
int lsf_history_alloc(lsf_history *h, int rx, int ry, int nchan);
int lsf_history_alloc_rej(lsf_history *h, int first_gen);
int lsf_history_grow(lsf_history *h, const int pad[4]);
void lsf_history_free(lsf_history *h);
void lsf_file_free(lsf_file *f);
int lsf_open_and_validate(const char *path, lsf_file *f, gboolean tolerate_edited_stack);
gchar *lsf_hash_float_data(const float *data, size_t n);
gchar *lsf_now_utc_iso(void);
int lsf_write_file(const char *path, lsf_file *old, const lsf_settings *s, lsf_history *h,
		fits *primary, const char *uid, int gen, GPtrArray *new_rows, int batch, gboolean backup,
		lsf_compression comp);
const char *lsf_compression_name(lsf_compression comp);
int lsf_export_master(const char *lsf_path, const char *out_path);
void lsf_strip_lsf_keys(fits *fit);

/* lsf_stack.c */
struct lsf_batch_args {
	sequence *seq;			// registered batch, on the grid
	sequence *src_seq;		// unregistered source (append), or NULL
	int *src_index;			// seq index -> src_seq index, or NULL
	const Homography *xfrm_base;	// init: maps the registration of src_seq onto the grid, or NULL
	const char *framing;		// framing of the registration ("max", "current"), or NULL
	lsf_settings *settings;
	lsf_history *hist;		// in: history (may be empty), out: merged
	gboolean init;			// first batch: freeze the constants
	int ref_image;			// reference frame index in seq (init)
	int batch;			// batch number
	int gen;			// generation of this commit
	long prev_frames;		// frames already in the ledger
	const char *seq_name;
	GPtrArray *rows;		// out: lsf_ledger_row* of the merged frames
	fits *primary;			// in/out: header metadata for HDU 0 (built from the reference frame on init)
	guint64 rejected[3][2], hist_rejected[3][2];
};
int lsf_process_batch(struct lsf_batch_args *args);
void lsf_ledger_row_free(gpointer p);
int lsf_build_primary(fits *primary, lsf_settings *s, const lsf_history *h);

#endif /* SRC_LSF_LSF_H_ */
