#ifndef SU2OPS_H
#define SU2OPS_H

#include <lattice.h>

#define ELM2x2(a,b)    (a) * 2 + (b)  //  used to get the element of matrices

typedef unsigned short SU2_color_idx;

typedef struct {
    work_data_type m[2 * 2];
} mtrx_2x2;

typedef struct {
    double m[4];
} mtrx_2x2_ck;


// void print_mtrx_2x2(const mtrx_2x2 * restrict u, const char *name, const unsigned short decimal_places);

// void copy_2x2(const mtrx_2x2_ck * restrict u, mtrx_2x2_ck* restrict u_copy);

void convert_from_ck(const mtrx_2x2_ck * restrict u_ck, mtrx_2x2 * restrict u);

// void set_null_2x2(mtrx_2x2_ck * u), set_identity_2x2(mtrx_2x2_ck * u);

// void accumulate_2x2(const mtrx_2x2_ck * u, mtrx_2x2_ck* acc);

// void subtraction_2x2(const mtrx_2x2_ck* u, const mtrx_2x2_ck* v, mtrx_2x2_ck* u_minus_v);

// double SU2_trace(const mtrx_2x2_ck *u),
//        determinant_2x2(const mtrx_2x2_ck *u);

// void herm_conj_2x2(const mtrx_2x2_ck* u, mtrx_2x2_ck* u_dagger);

// void mult_scalar_2x2(const mtrx_2x2_ck* restrict u, const double alpha, mtrx_2x2_ck* restrict alpha_times_u);

// void product_2x2(const mtrx_2x2_ck* u, const mtrx_2x2_ck* v, mtrx_2x2_ck* uv),
//      product_three_2x2(const mtrx_2x2_ck* u, const mtrx_2x2_ck* v, const mtrx_2x2_ck* w, mtrx_2x2_ck* uvw),
//      product_four_2x2(const mtrx_2x2_ck* u, const mtrx_2x2_ck* v, const mtrx_2x2_ck* w, const mtrx_2x2_ck* x, mtrx_2x2_ck* uvwx);

short SU2_projection(mtrx_2x2_ck* restrict u);

#endif