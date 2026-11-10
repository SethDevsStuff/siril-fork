#ifndef SRC_LSF_LSF_H_
#define SRC_LSF_LSF_H_

/* Living Stack Format (LSF), FITS profile v0.4: the constants and the types
 * shared by the LSF modules. "Section" in the comments refers to that
 * specification. */

#include <glib.h>

#include "core/settings.h"
#include "core/siril.h"
#include "registration/registration.h"
#include "stacking/stacking.h"

#define LSF_VERSION "0.4"
#define LSF_EXT_ACC "LSF_ACC"
#define LSF_EXT_COUNT "LSF_COUNT"
#define LSF_EXT_REJ "LSF_REJ"
#define LSF_EXT_SETTINGS "LSF_SETTINGS"
#define LSF_EXT_LEDGER "LSF_LEDGER"
#define LSF_EXT_WREF "LSF_WREF"
#define LSF_COUNT_MAX 65535
#define LSF_LOW_STACK_WARN 25
#define LSF_NEFF_MIN (1.0 + 1e-4)
/* the history test runs only where the history holds this many effective
 * samples (Section 10.6): as Siril's rejection, which stops at 4 */
#define LSF_HIST_MIN_NEFF 5.0
#define LSF_EXPAND_STEP 64	// default step of a grid expansion, px per side (Section 10.5)
#define LSF_FOCAL_TOLERANCE 0.05	// relative FOCALLEN difference still taken as the same optics (Section 7)
/* default least wfwhm or nbstars weight before normalization (weight.floor,
 * Section 10.3), 1 being the best frame of the first batch: a frame at or
 * beyond its worst still counts */
#define LSF_WEIGHT_FLOOR 0.01
/* a frame covering fewer reference stars (LSF_WREF) has an unreliable wFWHM,
 * which happens to mosaic panels away from the first batch */
#define LSF_WREF_MIN_COVER 10
/* least overlap with the Stack for -overlap_norm: pixels, and share of the
 * frame's data */
#define LSF_OVERLAP_MIN_PIX 10000
#define LSF_OVERLAP_MIN_FRAC 0.01
/* bins of -overlap_norm on initialization, px: the frames of the batch are
 * matched on the bins they share */
#define LSF_OVERLAP_BIN 16
/* a frame under max framing keeps at least this share of its area on the grid */
#define LSF_MIN_LANDED 0.9
#define LSF_FEATHER_MAX 2000	// as the stack command's -feather=, in pixels of its masks (1/10 scale)

/* how the frozen wfwhm and nbstars constants scale the frames (weight.model,
 * Section 10.3): RANGE, Siril's formula between the first batch's best and
 * worst (its span multiplied by weight.range), a frame beyond them clamped;
 * INVERSE, unclamped: (best wFWHM / wFWHM)^2, Siril's formula with an
 * infinite worst, or (stars / most stars)^2. Files start with RANGE, and an
 * append may switch them to INVERSE. */
enum { LSF_WMODEL_RANGE, LSF_WMODEL_INVERSE };

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
	/* exposure map (Section 3), the EXP planes of LSF_ACC, NULL when the file
	 * has none: sum of EXPTIME over the accepted samples, seconds */
	float *exp[3];
	GSList *spill;			// the planes in temporary files (lsf_spill), when they do not fit in memory (Section 10.10)
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
	int oversample;			// grid.oversample, the drizzle scale, 1 without drizzle
	gboolean drizzled;		// the frames are drizzled by LSF, each sample weighted (Section 10.3): weight.driz_norm is set, and cfa.mode is CFADRIZ on a CFA sensor
	gboolean undistorted;		// grid.undistorted: every frame is undistorted before it is merged
	/* registration */
	transformation_type transform;	// align.transform
	opencv_interpolation interpolation;	// align.interpolation
	int min_pairs;			// align.min_pairs
	gchar *findstar_json;		// align.findstar
	disto_params disto;		// align.disto, the distortion model of the appends (undistorted grid only)
	/* rejection */
	rejection rej;			// reject.method
	double sig[2];			// reject.params
	double hist_low, hist_high, hist_margin;
	/* weighting, with constants frozen from the first batch (§10.3) */
	weightingType weighting;	// weight.scheme
	double noise_norm[3];		// per channel
	double fwhm_min, fwhm_max, fwhm_norm;
	int weight_model;		// weight.model, LSF_WMODEL_* (wfwhm and nbstars)
	double weight_floor;		// weight.floor, the least weight before the norm
	double weight_range;		// weight.range, multiplier of the first batch's span (1: Siril's)
	double stars_min, stars_max, stars_norm;
	double driz_norm[3];		// weight.driz_norm, Siril's max_weight of the first batch, per channel (drizzled files)
	gboolean driz_flat;		// weight.driz_flat: the first batch's drizzle weights, and so driz_norm, include a master flat (Section 13: frames drizzled by Siril's registration)
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

/* LSF_WREF (Section 10.3): the stars of the reference frame of the first
 * batch, on the grid, frozen. Every frame's wFWHM is measured against them,
 * so that it does not depend on the batch, the registration or the depth of
 * the Stack. FITS pixel coordinates of the current grid. */
typedef struct {
	int n;
	double *x, *y;
} lsf_wref;

/* values from the existing file needed to rewrite it */
typedef struct {
	gchar *path;
	gchar uid[64];
	int gen;
	int nchan, rx, ry;
	lsf_settings settings;
	lsf_history hist;
	lsf_wref wref;			// LSF_WREF, n = 0 when the file has none
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
	gboolean has_exp;		// LSF_ACC holds the EXP planes of the exposure map
	/* EXPTIME of every ledger row when they all have the same, known one, so
	 * that the exposure map is COUNT times it; NAN otherwise */
	double ledger_exptime;
	/* DRZ_PIXF and DRZ_KERN of the last ledger row, the defaults of an append
	 * to a drizzled file (Section 6): NAN and "" when the ledger has none */
	double ledger_pixfrac;
	gchar ledger_kernel[16];
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
	gchar *reg_dist;	// REG_DIST: the SIP cards the frame was undistorted with, or NULL
	double norm_sc[3], norm_off[3];
	gchar norm_mth[3][13];
	int feather;		// FEATHER: the feathering distance the frame was merged with, px, 0 for none
	float pixfrac;		// DRZ_PIXF: the drizzle pixel fraction, drizzled files only
	gchar kernel[9];	// DRZ_KERN: the drizzle kernel, empty when the frame was not drizzled
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

#endif /* SRC_LSF_LSF_H_ */
