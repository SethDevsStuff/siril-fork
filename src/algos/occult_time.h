#ifndef SRC_ALGOS_OCCULT_TIME_H_
#define SRC_ALGOS_OCCULT_TIME_H_

struct occultation_args {
	int start_ind_inseq;	// Index of the first image of the pulse in the original sequence
	int start_ind;	// Index of the first image of the pulse in the sorted list
	int pls_nbr;	// Number of usefull images in the 100ms pulse
	double sum_flux;	// Total flux during the pulse
	double delay_comp;	// delay time between the real PPS and the forseen PPS
};

struct occ_res {
	double median_seq;	// median over the sequence
	double std_seq;	// Std over the sequence
	double exposure;	// Computed exposure time
	int th_pls_nbr;	// Theorical pulse number in the sequence
	int det_pulses;	// Number of detected pulses
	double hi_val;	// Maximum value in the sequence
	double lo_val;	// Minimum value in the sequence
	int valid_images;	// Number of valid images in the sequence
};

gpointer occultation_worker(gpointer arg);

int time_resync(struct light_curve_args *lcargs);

#endif /* SRC_ALGOS_SRC_ALGOS_OCCULT_TIME_H__H_ */