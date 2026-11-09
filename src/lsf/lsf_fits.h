#ifndef SRC_LSF_LSF_FITS_H_
#define SRC_LSF_LSF_FITS_H_

/* LSF file I/O: validation (Section 8), reading, writing and the commit
 * protocol (Sections 10.2, 10.10) */

#include <stddef.h>

#include <glib.h>

#include "core/siril.h"
#include "lsf/lsf.h"

/* what a user interface shows of an LSF file, read by lsf_read_info without
 * the planes or the checksums, so fast enough to run when a file is chosen */
typedef struct {
	lsf_file file;			// as lsf_open_and_validate, without hist and stack_data
	int failed_check;		// 0, or the check of Section 8 that failed (1 to 4, 6, 7)
	gchar *problem;			// what failed, NULL when the checks pass
	gchar *no_append;		// why this version cannot append to the file, or NULL
	double exposure;		// sum of EXPTIME over the ledger, s
	gchar *first_date, *last_date;	// DATE_OBS range of the ledger, NULL when unknown
	gchar *last_added;		// ADDED_UTC of the last ledger row, NULL when unknown
	gint64 file_bytes;
} lsf_info;

/* one ledger row as the LSF Manager shows it (Section 9) */
typedef struct {
	int batch, gen;
	gchar *filename, *date_obs, *filter;
	double exptime;
	double weight[3];		// per channel
	float fwhm, wfwhm, noise[3];
	int nstars;
} lsf_ledger_view;

/* results of lsf_write_file: the two replace failures leave a file the
 * user has to deal with (Section 10.10) */
enum { LSF_WRITE_OK, LSF_WRITE_FAILED, LSF_WRITE_KEPT_TMP, LSF_WRITE_LEFT_AS_BAK };

int lsf_history_alloc(lsf_history *h, int rx, int ry, int nchan);
int lsf_history_alloc_rej(lsf_history *h, int first_gen);
int lsf_history_alloc_exp(lsf_history *h);
guint64 lsf_history_bytes(guint64 npix, gboolean rej, gboolean exp);
gboolean lsf_accumulator_spills(guint64 npix);
int lsf_history_grow(lsf_history *h, const int pad[4]);
void lsf_history_free(lsf_history *h);
void lsf_file_free(lsf_file *f);
void lsf_wref_free(lsf_wref *w);
long lsf_read_ledger_weights(const lsf_file *f, float **wfwhm, int **nstars, double **weight);
int lsf_open_and_validate(const char *path, lsf_file *f, gboolean tolerate_edited_stack);
int lsf_open_and_validate_msg(const char *path, lsf_file *f, gboolean tolerate_edited_stack, gchar **problem);
int lsf_read_info(const char *path, lsf_info *info);
void lsf_info_free(lsf_info *info);
GPtrArray *lsf_read_ledger_view(const lsf_info *info);
void lsf_ledger_view_free(gpointer p);
GPtrArray *lsf_read_history(const char *path);
gchar *lsf_hash_float_data(const float *data, size_t n);
gchar *lsf_now_utc_iso(void);
int lsf_write_file(const char *path, lsf_file *old, const lsf_settings *s, lsf_history *h, const lsf_wref *wref,
		fits *primary, const char *uid, int gen, GPtrArray *new_rows, int batch, gboolean backup,
		lsf_compression comp, gchar **replace_msg);
const char *lsf_compression_name(lsf_compression comp);
int lsf_export_master(const char *lsf_path, const char *out_path);
void lsf_strip_lsf_keys(fits *fit);
void lsf_read_view_uid(const fits *fit, char *uid, size_t n);

#endif /* SRC_LSF_LSF_FITS_H_ */
