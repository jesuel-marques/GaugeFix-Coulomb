#ifndef SU2OPS_H
#define SU2OPS_H

#include <types.h>


// 2x2 types

typedef unsigned short MtrxIdx2;            //  2x2 matrix index

typedef struct {
    Scalar m[2 * 2];
} Mtrx2x2;                                  /*  Struct for 2x2 matrices */

typedef struct {
    double m[4];
} Mtrx2x2CK;                                /*  Struct for 2x2 matrices in the Cayley-
                                                Klein form 
                                                u = m[0] SU2_identity 
                                                    + i sum_i=1^3 m[i]sigma[i] */


#define ELM_2X2(a,b)    (a) * 2 + (b)  //  used to get the element of matrices

#define LOOP_2_CK(a)    for(a = 0; a < 4; a++)
#define LOOP_2_CK_i(a)  for(a = 1; a < 4; a++)

#define LOOP_2X2(a, b)  for(a = 0; a < 3; a++)  \
                        for(b = 0; b < 3; b++)

//void printMtrx2x2(const Mtrx2x2 * restrict u, const char *name, const unsigned short decimal_places);

void copy2x2(const Mtrx2x2CK * restrict u, 
                    Mtrx2x2CK* restrict u_copy);

void convertFromCK(const Mtrx2x2CK * restrict u_ck, 
                           Mtrx2x2 * restrict u);

//void setNull2x2(Mtrx2x2CK * u), setIdentity2x2(Mtrx2x2CK * u);

//void accumulate2x2(const Mtrx2x2CK * u, 
//                          Mtrx2x2CK* acc);
//void subtraction2x2(const Mtrx2x2CK* u, 
//                     const Mtrx2x2CK* v, 
//                           Mtrx2x2CK* u_minus_v);

//Scalar SU2Trace(const Mtrx2x2CK *u);

Scalar determinant2x2(const Mtrx2x2CK *u);
//void hermConj2x2(const Mtrx2x2CK* u, 
//                         Mtrx2x2CK* u_dagger);

void multScalar2x2(const Mtrx2x2CK * restrict u, 
                     const Scalar alpha, 
                           Mtrx2x2CK * restrict alpha_times_u);

//void product2x2(const Mtrx2x2CK* u, 
//                 const Mtrx2x2CK* v, 
//                       Mtrx2x2CK* uv),
//productThree2x2(const Mtrx2x2CK* u, 
//                  const Mtrx2x2CK* v, 
//                  const Mtrx2x2CK* w, 
//                        Mtrx2x2CK* uvw),
//productFour2x2(const Mtrx2x2CK* u, 
//                 const Mtrx2x2CK* v,
//                 const Mtrx2x2CK* w, 
//                 const Mtrx2x2CK* x, 
//                       Mtrx2x2CK* uvwx);

short projectSU2(Mtrx2x2CK* restrict u);

#endif  //SU2OPS_H