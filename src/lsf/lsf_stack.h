#ifndef SRC_LSF_LSF_STACK_H_
#define SRC_LSF_LSF_STACK_H_

/* LSF merge of a batch into the accumulator (Sections 10.2 to 10.6) */

#include <glib.h>

#include "core/siril.h"
#include "lsf/lsf.h"

struct lsf_batch_args {
	sequence *seq;			// registered batch, on the grid
	sequence *src_seq;		// unregistered source (append), or NULL
	int *src_index;			// seq index -> src_seq index, or NULL
	const Homography *xfrm_base;	// init: maps the registration of src_seq onto the grid, or NULL
	const char *framing;		// framing of the registration ("max", "current"), or NULL
	const char *reg_method;		// append: "stars" or "wcs" (Section 10.4), NULL for stars
	const disto_params *disto;	// distortion model the frames were undistorted with, or NULL
	lsf_settings *settings;
	lsf_history *hist;		// in: history (may be empty), out: merged
	gboolean init;			// first batch: freeze the constants
	int ref_image;			// reference frame index in seq (init)
	int batch;			// batch number
	int gen;			// generation of this commit
	long prev_frames;		// frames already in the ledger
	double prev_exptime;		// append: lsf_file.ledger_exptime, for the exposure map
	const char *seq_name;
	GPtrArray *rows;		// out: lsf_ledger_row* of the merged frames
	fits *primary;			// in/out: header metadata for HDU 0 (built from the reference frame on init)
	lsf_wref *wref;			// LSF_WREF (Section 10.3), on the current grid: filled on init
	gboolean keep_seq_file;		// the selection of seq was changed in memory: do not write its .seq file
	int max_threads;		// thread budget of the merge, at least 1
	gboolean overlap_norm;		// append: normalize against the Stack on the overlap
	int feather;			// feathering distance of the frames' edges, px, 0 for none
	guint64 rejected[3][2], hist_rejected[3][2];
};

int lsf_process_batch(struct lsf_batch_args *args);
int lsf_solve_levels(int n, int ref, size_t **Nij, double **Mij, gboolean additive, double *coeffs, gboolean *reached);
void lsf_ledger_row_free(gpointer p);
const char *lsf_frame_mismatch(const fits *primary, const fits *frame);
int lsf_build_primary(fits *primary, lsf_settings *s, const lsf_history *h);

#endif /* SRC_LSF_LSF_STACK_H_ */
