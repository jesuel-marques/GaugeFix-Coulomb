#ifndef SU3OPS_H
#define SU3OPS_H

#include <types.h>

#define ELEM_3X3(a, b)    (a) * Nc + (b)  //  used to get the element of matrices

#define LOOP_3(a)         for ( a = 0; a < Nc; a++)

#define LOOP_3X3(a, b)    for ( a = 0; a < Nc; a++)  \
                                for ( b = 0; b < Nc; b++)


void print_matrix_3x3(const Mtrx3x3 *u, 
                      const char *name, 
                      const unsigned short decimal_places);

void copy_3x3(const Mtrx3x3 *u, 
                    Mtrx3x3 *u_copy);

void set_null_3x3    (Mtrx3x3 *u),
     set_identity_3x3(Mtrx3x3 *u);

void accumulate_3x3(const Mtrx3x3 *u, 
                          Mtrx3x3 *acc);

void subtraction_3x3(const Mtrx3x3 *u, 
                     const Mtrx3x3 *v, 
                           Mtrx3x3 *u_minus_v);

Scalar trace_3x3      (const Mtrx3x3 *u),
       determinant_3x3(const Mtrx3x3 *u);

void herm_conj_3x3(const Mtrx3x3 *u, 
                         Mtrx3x3 *u_dagger);

void subtraction_herm_conj_trless_3x3(const Mtrx3x3 * restrict u, 
                                               Mtrx3x3 * restrict u_minus_udagger);

void mult_by_scalar_3x3(const Scalar alpha, 
                        const Mtrx3x3 *u, 
                              Mtrx3x3 *alpha_times_u);
void subst_mult_scalar_3x3(const Scalar alpha, 
                                 Mtrx3x3 *u);

void prod_3x3(const Mtrx3x3 *u, const Mtrx3x3 *v, Mtrx3x3 *uv);

void prod_three_3x3(const Mtrx3x3 *u, 
                    const Mtrx3x3 *v, 
                    const Mtrx3x3 *w, 
                          Mtrx3x3 *uvw),
     prod_four_3x3 (const Mtrx3x3 *u, 
                    const Mtrx3x3 *v, 
                    const Mtrx3x3 *w, 
                    const Mtrx3x3 *x, 
                          Mtrx3x3 *uvwx);

void accum_left_prod_3x3 (const Mtrx3x3 *g, 
                                Mtrx3x3 *acc_prod),
     accum_right_prod_3x3(      Mtrx3x3 *acc_prod, 
                          const Mtrx3x3 *g);

void accum_prod_SU2_3x3(const Mtrx2x2CK * restrict x_ck, 
                              Mtrx3x3 * restrict g, 
                              MtrxIdx3 a, 
                              MtrxIdx3 b);

void prod_vuwdagger_3x3(const Mtrx3x3 * restrict v, 
                        const Mtrx3x3 * restrict u, 
                        const Mtrx3x3 * restrict w, 
                              Mtrx3x3 * restrict vuwdagger );

short power_3x3_binomial(Mtrx3x3 * restrict A, 
                         const double omega, 
                         Mtrx3x3 * restrict A_to_omega);


double inverse_3x3(const Mtrx3x3 * restrict a, 
                         Mtrx3x3 * restrict a_inv);

short projection_SU3(Mtrx3x3 *x);


void decompose_algebra_SU3(const Mtrx3x3 *a, 
                                 MtrxSU3Alg *a_components);

void SU3_CabbiboMarinari_projection(Mtrx3x3 * restrict w, 
                                    Mtrx3x3 * restrict total_update);

#endif      //SU3OPS_H