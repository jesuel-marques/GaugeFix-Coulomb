#ifndef SU3OPS_H
#define SU3OPS_H

#include <types.h>
#include <SU2_ops.h>

#define Nc  3                   //  Number of colors

// 3x3 types

typedef unsigned short MtrxIdx3;            //  3x3 matrix index
typedef unsigned short SU3AlgIdx;           //  SU(3) algebra index

typedef struct {
   complex double m[Nc * Nc];
} Mtrx3x3Dbl;                               /*  Struct for 3x3 matrices in double   */

typedef struct {
    complex float m[Nc * Nc];
} Mtrx3x3Flt;                               /*  Struct for 3x3 matrices in float    */

typedef struct {
   Scalar m[Nc];
} Vec3;                                     /*  Struct for 3 component vectors in color
                                                space. */

typedef struct {
    Scalar m[(Nc * Nc - 1) + 1];
} MtrxSU3Alg;                               /*  Struct for SU(3) algebra components. */

typedef Mtrx3x3Dbl WorkMtrx;                //  Calculations internally done in double
typedef WorkMtrx Mtrx3x3;                   //  Default matrix is matrix of type Work

typedef enum {R, S, T} Submtrx;             /*  Submtrx provides an alias for the 
                                                indices of Cabbibo-Marinari subdivision 
                                                of an SU(3) matrix.

/*
This is a way of embedding SU(2) matrices in an SU(3) matrix. For each of the indices,
the X marks where the SU(2) matrix goes and the rest is kept as either a 0 or a 1, 
as shown.

R ->  X X 0
      X X 0
      0 0 1

S ->  X 0 X
      0 1 0
      X 0 X

T ->  1 0 0
      0 X X
      0 X X

This is used in this program for the updates to the local gauge transformations. The 
updates are performed for each submatrix in turn.
*/


#define ELEM_3X3(a, b)    (a) * Nc + (b)  //  used to get the element of matrices

#define LOOP_3(a)         for( a = 0; a < Nc; a++)

#define LOOP_3X3(a, b)    for( a = 0; a < Nc; a++)  \
                              for( b = 0; b < Nc; b++)


void printMatrix3x3(const Mtrx3x3 *u);

void copy3x3(const Mtrx3x3 *u, 
             Mtrx3x3 *u_copy);

void setNull3x3    (Mtrx3x3 *u),
     setIdentity3x3(Mtrx3x3 *u);

void accumulate3x3(const Mtrx3x3 *u, 
                   Mtrx3x3 *acc);

void subtraction3x3(const Mtrx3x3 *u, 
                    const Mtrx3x3 *v, 
                    Mtrx3x3 *u_minus_v);

Scalar trace3x3      (const Mtrx3x3 *u),
       determinant3x3(const Mtrx3x3 *u);

void hermConj3x3(const Mtrx3x3 *u, 
                 Mtrx3x3 *u_dagger);

void subtractionHermConjTrless3x3(const Mtrx3x3 * restrict u, 
                                  Mtrx3x3 * restrict u_minus_udagger);

void multByScalar3x3(const Scalar alpha, 
                     const Mtrx3x3 *u, 
                     Mtrx3x3 *alpha_times_u);
void substMultScalar3x3(const Scalar alpha, 
                        Mtrx3x3 *u);

void prod3x3(const Mtrx3x3 *u, const Mtrx3x3 *v, Mtrx3x3 *uv);

void prodThree3x3(const Mtrx3x3 *u, 
                  const Mtrx3x3 *v, 
                  const Mtrx3x3 *w, 
                  Mtrx3x3 *uvw),
     prodFour3x3 (const Mtrx3x3 *u, 
                  const Mtrx3x3 *v, 
                  const Mtrx3x3 *w, 
                  const Mtrx3x3 *x, 
                  Mtrx3x3 *uvwx);

void accumLeftProd3x3 (const Mtrx3x3 *g, 
                       Mtrx3x3 *acc_prod),
     accumRightProd3x3(Mtrx3x3 *acc_prod, 
                       const Mtrx3x3 *g);

void accumProdSU2_3x3(const Mtrx2x2CK * restrict x_ck, 
                      Mtrx3x3   * restrict g, 
                      MtrxIdx3 a, 
                      MtrxIdx3 b);

void prod_vuwdagger3x3(const Mtrx3x3 * restrict v, 
                       const Mtrx3x3 * restrict u, 
                       const Mtrx3x3 * restrict w, 
                       Mtrx3x3 * restrict vuwdagger );

short power3x3Binomial(const Mtrx3x3 * restrict A, 
                       const double omega, 
                       Mtrx3x3 * restrict A_to_omega);

void updateSubLosAlamos(Mtrx3x3 * restrict w, Submtrx sub);

double inverse3x3(const Mtrx3x3 * restrict a, 
                  Mtrx3x3 * restrict a_inv);

short projectSU3(Mtrx3x3 *x);

void decomposeAlgebraSU3(const Mtrx3x3 *a, 
                         MtrxSU3Alg *a_components);

void projectSU3CabbiboMarinari(Mtrx3x3 * restrict w, 
                               Mtrx3x3 * restrict total_update);

#endif      //SU3OPS_H