#include "io/sequence.h"
#include "registration/lena.h"
#include <stdlib.h>
#include <sys/time.h>

float thr_pc = 0.7f;

void lena_set_thr_pc(float pc) {
    thr_pc = pc;
}

void lena_fingerprint_hz(fits *img, WORD max, int layer, int x0, int y0, int w, int h, lena_fingerprint *fp) {
    fp->type = img->type;
    WORD thr = (WORD) max*thr_pc;

    if (fp->type == DATA_USHORT) {
        if (fp->pdata) free(fp->pdata);
        fp->pdata = (WORD*) calloc(h, sizeof(WORD));
        fp->len = h;
        for (int y = y0; y < y0+h; y++) {
            for (int x = x0; x < x0 +w; x++) {
                WORD v = img->pdata[layer][x + y * img->rx];
                if (v>=thr)
                    fp->pdata[y - y0] += v;
            }
        }
    }
}

void lena_fingerprint_vt(fits *img, WORD max, int layer, int x0, int y0, int w, int h, lena_fingerprint *fp) {
    fp->type = img->type;
    WORD thr = (WORD) max*thr_pc;

    if (fp->type == DATA_USHORT) {
        if (fp->pdata) free(fp->pdata);
        fp->pdata = (WORD*) calloc(w, sizeof(WORD));
        fp->len = w;
        for (int y = y0; y < y0+h; y++) {
            for (int x = x0; x < x0 +w; x++) {
                WORD v = img->pdata[layer][x + y * img->rx];
                if (v>=thr)
                    fp->pdata[x - x0] += v;
            }
        }
    }
}

void lena_fingerprint_set_max(fits *img, int layer, lena_fingerprint *hz, lena_fingerprint *vt)
{
    if (hz->type == DATA_USHORT) {
        WORD max = 0;
        
        for (int y = 0; y < img->ry; y++) {
            for (int x = 0; x < img->rx; x++) {
                WORD v = img->pdata[layer][x + y * img->rx];
                if (v>max) max = v;
            }
        }

        hz->pdata_max = vt->pdata_max = max;
    }
}

void lena_fingerprint_set_thr(fits *img, int layer, lena_fingerprint *hz, lena_fingerprint *vt)
{
    if (hz->type == DATA_USHORT) {
        WORD thr_hz = 0, thr_vt = 0;

        for(int i = 0; i < hz->len; i++)
            thr_hz += hz->pdata[i];

        for(int i = 0; i < vt->len; i++)
            thr_vt += vt->pdata[i];

        hz->pdata_thr_hz = thr_hz;
        vt->pdata_thr_vt = thr_vt;

        vt->pdata_thr_hz = thr_hz;
        hz->pdata_thr_vt = thr_vt;
    }
}


float lena_diff_fingerprint(lena_fingerprint *f0, lena_fingerprint *f1) {
    float diff = 0;
    int r;

    if ((f0->type == DATA_USHORT) && (f1->type == DATA_USHORT))
        for(int i = 0; i < f0->len; i++) {
            r = f0->pdata[i] - f1->pdata[i];
    
            if (r<0) r = -r;
            diff += r;
        }

    return(diff);
}

void lena_min_fingerprint_vt(lena_fingerprint *vt_a, int layer, fits *l_b, int dx, int x0, int y0, int w, int h,
                             int *min_x, float *min_dt) {
    float mdt = -1;
    int mx = -1;
    lena_fingerprint vt_b = { 0 };
        
    WORD max = vt_a->pdata_max;
    int f_min_x = x0 - dx;
    int f_max_x = x0 + dx;

    if (f_min_x<0) f_min_x = 0;
    if (f_max_x>= l_b->rx) f_max_x = l_b->rx -1;

    for (int rx = f_min_x; rx < f_max_x; rx++) {
        lena_fingerprint_vt(l_b, max, layer, rx, y0, w, h, &vt_b);
        float dt = lena_diff_fingerprint(vt_a, &vt_b);
              
        if ((mdt>dt) || (mdt == -1)) {
            mdt = dt;
            mx = rx - x0;
        }
    }

    *min_x = mx;
    *min_dt = mdt;
}

void lena_min_fingerprint_hz(lena_fingerprint *hz_a, int layer, fits *l_b, int dy, int x0, int y0, int w, int h,
                             int *min_y, float *min_dt) {
    int my = -1;
    float mdt = -1;
    lena_fingerprint hz_b = { 0 };
        
    WORD max = hz_a->pdata_max;
    int f_min_y = y0 - dy;
    int f_max_y = y0 + dy;

    if (f_min_y<0) f_min_y = 0;
    if (f_max_y>= l_b->ry) f_max_y = l_b->ry - 1;


    for (int ry = f_min_y; ry < f_max_y; ry++) {
        lena_fingerprint_hz(l_b, max, layer, x0, ry, w, h, &hz_b);
        float dt = lena_diff_fingerprint(hz_a, &hz_b);

        if ((mdt>dt) || (mdt == -1)) {
            mdt = dt;
            my = ry - y0;
        }
    }

    *min_y = my;
    *min_dt = mdt;
}

void lena_register_one(reg_lena *args) {
    lena_fingerprint *vt_a = args->vt_a;
    lena_fingerprint *hz_a = args->hz_a;
    fits *l_b = args->l_b;
    int layer = args->layer;
    int scope_x = args->scope_x;
    int scope_y = args->scope_y;
    float *best_dt = &args->best_dt;
    int *best_dx = &args->best_dx;
    int *best_dy = &args->best_dy;
    int x0 = args->x0;
    int y0 = args->y0;
    int w  = args->w;
    int h  = args->h;
    
    float dt_x = 0, dt_y = 0;
    int dx = 0, dy = 0;

    *best_dt = *best_dx = *best_dy = -1;

    int turns = 2;

    while (turns--) {
        lena_min_fingerprint_vt(vt_a, layer, l_b, scope_x, x0, y0 + dy, w, h, &dx, &dt_x);
        lena_min_fingerprint_hz(hz_a, layer, l_b, scope_y, x0 + dx, y0, w, h, &dy, &dt_y);
    }
    *best_dt = dt_x + dt_y;
    *best_dx = dx;
    *best_dy = dy;

    args->valid = 1; // (score<=0.75);
}