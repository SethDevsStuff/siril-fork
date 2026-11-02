#ifndef SRC_LSF_LSF_COMMANDS_H_
#define SRC_LSF_LSF_COMMANDS_H_

/* LSF commands and the same jobs started by the GUI. The process_lsf_*
 * command functions are declared with the other commands, in core/command.h */

#include <glib.h>

#include "core/siril.h"
#include "lsf/lsf.h"
#include "lsf/lsf_fits.h"
#include "registration/registration.h"

typedef enum { LSF_CMD_INIT, LSF_CMD_APPEND, LSF_CMD_EXPORT, LSF_CMD_UTIL, LSF_CMD_MAPS } lsf_cmd_mode;

/* maps of lsf_maps (Section 10.8) */
enum { LSF_MAP_NEFF = 1 << 0, LSF_MAP_NOISE = 1 << 1, LSF_MAP_SIGMA = 1 << 2, LSF_MAP_SNR = 1 << 3,
	LSF_MAP_REJLO = 1 << 4, LSF_MAP_REJHI = 1 << 5, LSF_MAP_EXP = 1 << 6 };

typedef struct lsf_job_result lsf_job_result;
/* called on the main thread when the job is over, after the processing
 * thread is released, so it may start another job */
typedef void (*lsf_end_func)(const lsf_job_result *res);

struct lsf_job_result {
	lsf_cmd_mode mode;
	int retval;		// 0 on success
	gchar *path;		// the LSF file
	gchar *out;		// util and append: the new file written instead of rewriting path, or NULL
	gchar *export_out;	// the plain image written with it, as named on disk, or NULL
	gchar *message;		// util: what was found or done, for the user, or NULL
	gboolean stack_modified;	// util: check 8 found the Stack edited outside Siril
	int write_status;	// LSF_WRITE_* of the file written; on a failed replace, message says what to do
	lsf_end_func end_func;
	gpointer end_data;	// of the caller
};

struct lsf_cmd_args {
	lsf_cmd_mode mode;
	gchar *seqname;
	gchar *path;
	gchar *out;		// export, the maps prefix, or for util and append the file to write instead of rewriting path
	lsf_settings settings;	// init: requested settings
	gboolean hist_low_set, hist_high_set, hist_margin_set;
	double hist_low, hist_high, hist_margin;
	gboolean noreg, backup, overwrite;
	gboolean rejmaps;		// init and append: -rejmaps (Section 5.4)
	framing_type framing;		// init and append: -framing=, FRAMING_MAX by default (Sections 10.1, 10.5)
	gboolean framing_set;
	gboolean expand_force;
	int expand_step;
	gboolean check, repair, adopt, info;	// util
	lsf_compression compression;	// -compress=, -nocompress, -decompress (Section 5.3)
	gboolean compression_set;	// append: keeps the compression of the file when not set
	guint maps;		// maps, LSF_MAP_* flags
	gchar *export_out;	// init and append: also export the new Stack to this plain image (GUI)
	gboolean load_result;	// init and append: the new Stack becomes the result image (GUI)
	gboolean rej_set;	// append: settings.rej and sig override the rejection stored in the file
	gboolean force_reg;	// append: -reg, register again even with a registration made against the file
	gboolean wcsreg;	// append: -wcsreg, register through the WCS of the frames and of the Stack (Section 10.4)
	gboolean weight_inverse;	// append: -weight_model=inverse, from this append on (Section 10.3)
	gboolean overlap_norm;		// append: -overlap_norm, normalization against the Stack on the overlap
	int feather;			// init and append: -feather=, px, 0 for none
	disto_params disto;	// append: -disto=, the distortion model of an undistorted grid; DISTO_UNDEF: the file's
	gchar *wcs_file;	// util: -wcs=, the plate-solved image or .wcs file whose WCS is imported (Section 8)
	int *image_indices;	// init and append: frames chosen by the GUI's stacking filters, or NULL
	int nb_images;
	int max_threads;	// thread budget of the job, com.max_thread by default
	gchar *message;		// out, util and a failed replace: passed to lsf_job_result
	gboolean stack_modified;	// out, util
	int write_status;	// out: LSF_WRITE_* of the file written
	lsf_end_func end_func;	// GUI: called when the job is over, or NULL
	gpointer end_data;
};

struct lsf_cmd_args *lsf_cmd_args_new(lsf_cmd_mode mode);
void lsf_cmd_args_free(struct lsf_cmd_args *a);
gchar *lsf_cmd_args_to_command(const struct lsf_cmd_args *a);
int lsf_start(struct lsf_cmd_args *a);
gchar *lsf_check_batch(const lsf_info *info, sequence *seq, gboolean *fatal);
int lsf_prepare_registration(struct registration_args *regargs, gboolean keep_type, gboolean keep_pairs,
		gboolean keep_findstar);
const gchar *lsf_registration_result(const struct registration_args *regargs);
gboolean lsf_binding_check(const sequence *seq, const lsf_file *f, const disto_params *model, int *layer,
		gchar **date, gchar **why);
gchar *lsf_append_preflight(const lsf_info *info, const char *wd, const char *out, int nb_frames, gboolean noreg,
		gboolean add_rejmaps, gboolean *fatal);
gboolean lsf_same_file(const char *p1, const char *p2);
gchar *lsf_init_preflight(sequence *seq, const char *wd, const char *path, int nb_frames, gboolean rejmaps,
		gboolean *fatal);

#endif /* SRC_LSF_LSF_COMMANDS_H_ */
