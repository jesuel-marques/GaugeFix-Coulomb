#ifndef SU2OPS_H
#define SU2OPS_H

#include <types.h>

#define ELM_2X2(a,b)    (a) * 2 + (b)  //  used to get the element of matrices

#define LOOP_2_CK(a)    for (a = 0; a < 4; a++)
#define LOOP_2_CK_i(a)  for (a = 1; a < 4; a++)

#define LOOP_2X2(a, b)  for ( a = 0; a < Nc; a++)  \
                        for ( b = 0; b < Nc; b++)

//void print_mtrx_2x2(const Mtrx2x2 * restrict u, const char *name, const unsigned short decimal_places);

void copy_2x2(const Mtrx2x2CK * restrict u, 
                    Mtrx2x2CK* restrict u_copy);

void convert_from_ck(const Mtrx2x2CK * restrict u_ck, 
                           Mtrx2x2 * restrict u);

//void set_null_2x2(Mtrx2x2CK * u), set_identity_2x2(Mtrx2x2CK * u);

//void accumulate_2x2(const Mtrx2x2CK * u, 
//                          Mtrx2x2CK* acc);

//void subtraction_2x2(const Mtrx2x2CK* u, 
//                     const Mtrx2x2CK* v, 
//                           Mtrx2x2CK* u_minus_v);

//Scalar SU2_trace(const Mtrx2x2CK *u);
Scalar determinant_2x2(const Mtrx2x2CK *u);

//void herm_conj_2x2(const Mtrx2x2CK* u, 
//                         Mtrx2x2CK* u_dagger);

void mult_scalar_2x2(const Mtrx2x2CK * restrict u, 
                     const Scalar alpha, 
                           Mtrx2x2CK * restrict alpha_times_u);

//void product_2x2(const Mtrx2x2CK* u, 
//                 const Mtrx2x2CK* v, 
//                       Mtrx2x2CK* uv),
//product_three_2x2(const Mtrx2x2CK* u, 
//                  const Mtrx2x2CK* v, 
//                  const Mtrx2x2CK* w, 
//                        Mtrx2x2CK* uvw),
//product_four_2x2(const Mtrx2x2CK* u, 
//                 const Mtrx2x2CK* v,
//                 const Mtrx2x2CK* w, 
//                 const Mtrx2x2CK* x, 
//                       Mtrx2x2CK* uvwx);

short SU2_projection(Mtrx2x2CK* restrict u);

#endif  //SU2OPS_H