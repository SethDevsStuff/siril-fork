#ifndef SRC_LSF_LSF_SPILL_H_
#define SRC_LSF_LSF_SPILL_H_

#include <stddef.h>

/* A zeroed block of memory backed by a temporary file, for an accumulator
 * larger than the memory Siril may use (Section 10.10): the system keeps in
 * memory the parts in use and writes the others to the file, which is removed
 * when the block is freed, or by the system if Siril stops. */
typedef struct lsf_spill lsf_spill;

lsf_spill *lsf_spill_new(const char *dir, size_t bytes);
void *lsf_spill_data(const lsf_spill *sp);
void lsf_spill_free(lsf_spill *sp);

#endif /* SRC_LSF_LSF_SPILL_H_ */
