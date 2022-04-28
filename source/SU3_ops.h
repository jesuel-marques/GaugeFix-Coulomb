#ifndef SU3OPS_H
#define SU3OPS_H

#include "lattice.h"

#define elm(a,b)    (a) * Nc + (b)  //  used to get the element of matrices

typedef unsigned short SU3_color_idx;
typedef unsigned short SU3_alg_idx;

typedef struct {
   work_data_type m[Nc];
} color_3_vec;

typedef struct {
    work_data_type m[Nc * Nc - 1 + 1];
} matrix_SU3_alg;

void print_matrix_3x3(const mtrx_3x3 *u, const char *name, const unsigned short decimal_places);

void copy_3x3(const mtrx_3x3 *u, mtrx_3x3 *u_copy);

void convert_in_work_3x3(const in_cfg_data_type *u_in, work_cfg_data_type *u_work),
     convert_work_out_3x3(const work_cfg_data_type *u_work, out_cfg_data_type *u_out);

void set_null_3x3    (mtrx_3x3 *u),
     set_identity_3x3(mtrx_3x3 *u);

void accumulate_3x3(const mtrx_3x3 *u, mtrx_3x3 *acc);

void subtraction_3x3(const mtrx_3x3 *u, 
                     const mtrx_3x3 *v, 
                                mtrx_3x3 *u_minus_v);

work_data_type trace_3x3(const mtrx_3x3 *u),
         determinant_3x3(const mtrx_3x3 *u);

void SU3_herm_conj(const mtrx_3x3 *u, mtrx_3x3 *u_dagger);

void mult_by_scalar_3x3(const work_data_type alpha, const mtrx_3x3 *u, mtrx_3x3 *alpha_times_u);
void subst_mult_scalar_3x3(const work_data_type alpha, mtrx_3x3 *u);

void prod_3x3(const mtrx_3x3 *u, const mtrx_3x3 *v, mtrx_3x3 *uv);

void prod_three_3x3(const mtrx_3x3 *u, const mtrx_3x3 *v, const mtrx_3x3 *w, mtrx_3x3 *uvw),
     prod_four_3x3(const mtrx_3x3 *u, const mtrx_3x3 *v, const mtrx_3x3 *w, const mtrx_3x3 *x, mtrx_3x3 *uvwx);

void accum_left_prod_3x3(const mtrx_3x3 *g, mtrx_3x3 *acc_prod),
     accum_right_prod_3x3(mtrx_3x3 *acc_prod, const mtrx_3x3 *g);

void projection_SU3(mtrx_3x3 *x);

void decompose_algebra_SU3(const mtrx_3x3 *a, matrix_SU3_alg *a_components);

#endif