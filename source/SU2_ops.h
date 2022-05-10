#ifndef SU2OPS_H
#define SU2OPS_H

#include "lattice.h"

typedef unsigned short SU2_color_idx;


typedef struct {
    double m[4];
} matrix_2x2_ck;

typedef struct {
    work_data_type m[2 * 2];
}matrix_2x2;

void print_matrix_2x2(const matrix_2x2 * restrict u, const char *name, const unsigned short decimal_places);

void copy_2x2(const matrix_2x2_ck * restrict u, matrix_2x2_ck* restrict u_copy);

void convert_from_ck(const matrix_2x2_ck * restrict u_ck, matrix_2x2 * restrict u);

// void set_null_2x2(matrix_2x2_ck * u), set_identity_2x2(matrix_2x2_ck * u);

// void accumulate_2x2(const matrix_2x2_ck * u, matrix_2x2_ck* acc);

// void subtraction_2x2(const matrix_2x2_ck* u, const matrix_2x2_ck* v, matrix_2x2_ck* u_minus_v);

// double SU2_trace(const matrix_2x2_ck *u),
//        determinant_2x2(const matrix_2x2_ck *u);

void herm_conj_2x2(const matrix_2x2_ck* u, matrix_2x2_ck* u_dagger);

void mult_scalar_2x2(const matrix_2x2_ck* restrict u, const double alpha, matrix_2x2_ck* restrict alpha_times_u);

// void product_2x2(const matrix_2x2_ck* u, const matrix_2x2_ck* v, matrix_2x2_ck* uv),
//      product_three_2x2(const matrix_2x2_ck* u, const matrix_2x2_ck* v, const matrix_2x2_ck* w, matrix_2x2_ck* uvw),
//      product_four_2x2(const matrix_2x2_ck* u, const matrix_2x2_ck* v, const matrix_2x2_ck* w, const matrix_2x2_ck* x, matrix_2x2_ck* uvwx);

void SU2_projection(matrix_2x2_ck* restrict u);

#endif