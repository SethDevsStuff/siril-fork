#ifndef H_LENA_H
#define H_LENA_H

#include "registration/registration.h"


typedef struct {
    data_type type;
    float     *fpdata;
    float     fpdata_thr;

    WORD      *pdata;
    WORD      pdata_thr_hz;
    WORD      pdata_thr_vt;
    WORD      pdata_max;
    
    int       len;
} lena_fingerprint;

typedef struct {
    lena_fingerprint *vt_a;
    lena_fingerprint *hz_a;
    fits *l_b;
    int layer;
    int x0, y0, w, h;
    int scope_x;
    int scope_y;
    float best_dt;
    int best_dx;
    int best_dy;
    int valid;
} reg_lena;

void lena_set_thr_pc(float pc);

void lena_fingerprint_hz(fits *img, WORD max, int layer, int x0, int y0, int w, int h, lena_fingerprint *fp);
void lena_fingerprint_vt(fits *img, WORD max, int layer, int x0, int y0, int w, int h, lena_fingerprint *fp);

void lena_register_one(reg_lena *args);

void lena_fingerprint_set_max(fits *img, int layer, lena_fingerprint *hz, lena_fingerprint *vt);
void lena_fingerprint_set_thr(fits *img, int layer, lena_fingerprint *hz, lena_fingerprint *vt);


#endif