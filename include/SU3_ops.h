/*
    Header to SU3_ops.c, which performs 3x3 matrix operations.

    Copyright (C) 2023  Jesuel Marques

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Contact: jesuel.leal@usp.br

 */

#ifndef SU3OPS_H
#define SU3OPS_H

#include <SU2_ops.h>
#include <types.h>

/* Number of colors */
#define Nc 3

// 3x3 types

/*  3x3 matrix index */
typedef unsigned short MtrxIdx3;

/* SU(3) algebra index */
typedef unsigned short SU3AlgIdx;

/* 3x3 complex matrices in double */
typedef struct Mtrx3x3Dbl {
    complex double m[Nc * Nc];
} Mtrx3x3Dbl;

/* 3x3 complex matrices in float */
typedef struct Mtrx3x3Flt {
    complex float m[Nc * Nc];
} Mtrx3x3Flt;

/* 3 component vectors in color space. */
typedef struct Vec3 {
    Scalar m[Nc];
} Vec3;

/* SU(3) algebra matrix with Nc²-1 components. */
typedef struct MtrxSU3Alg {
    Scalar m[(Nc * Nc - 1) + 1]; /* one component added to keep the numbering neat
                                    the zeroth component should not be used */
} MtrxSU3Alg;

/* The internal working precision of the calculations is double */
typedef Mtrx3x3Dbl Mtrx3x3Work;

/* The default matrix is a matrix of type Mtrx3x3Work */
typedef Mtrx3x3Work Mtrx3x3;

/*  Submtrx provides an alias for the indices of
    Cabbibo-Marinari subdivision of an SU(3) matrix. */
typedef enum { R,
               S,
               T } Submtrx;

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

typedef enum { ONE,
               TWO_PI_OVER_THREE,
               FOUR_PI_OVER_THREE } CenterElement;

/* used to get the element of 3x3 matrices */
#define ELEM_3X3(a, b) (a) * Nc + (b)

/* used to loop through color indices for 3x3 matrices or 3-vectors */
#define LOOP_3(a) for (a = 0; a < Nc; a++)

/* used to loop through color indices for 3x3 matrices */
#define LOOP_3X3(a, b)       \
    for (a = 0; a < Nc; a++) \
        for (b = 0; b < Nc; b++)

void printMatrix3x3(const Mtrx3x3 *u);

void copy3x3(const Mtrx3x3 *u, Mtrx3x3 *u_copy);

void setNull3x3(Mtrx3x3 *u),
    setIdentity3x3(Mtrx3x3 *u),
    setSU3Random(Mtrx3x3 *u);
void accumulate3x3(const Mtrx3x3 *u, Mtrx3x3 *acc);

void subtraction3x3(const Mtrx3x3 *u, const Mtrx3x3 *v, Mtrx3x3 *u_minus_v);

Scalar trace3x3(const Mtrx3x3 *u),
    determinant3x3(const Mtrx3x3 *u);

void hermConj3x3(const Mtrx3x3 *u, Mtrx3x3 *u_dagger);

void subtractionHermConjTrless3x3(const Mtrx3x3 *restrict u,
                                  Mtrx3x3 *restrict u_minus_udagger);

void multByScalar3x3(const Scalar alpha, const Mtrx3x3 *u, Mtrx3x3 *alpha_times_u);
void substMultScalar3x3(const Scalar alpha, Mtrx3x3 *u);

void prod3x3(const Mtrx3x3 *u, const Mtrx3x3 *v, Mtrx3x3 *uv);

void prodThree3x3(const Mtrx3x3 *u, const Mtrx3x3 *v, const Mtrx3x3 *w, Mtrx3x3 *uvw),
    prodFour3x3(const Mtrx3x3 *u,
                const Mtrx3x3 *v,
                const Mtrx3x3 *w,
                const Mtrx3x3 *x,
                Mtrx3x3 *uvwx);

void accumLeftProd3x3(const Mtrx3x3 *g, Mtrx3x3 *acc_prod),
    accumRightProd3x3(Mtrx3x3 *acc_prod, const Mtrx3x3 *g);

void accumProdSU2_3x3(const Mtrx2x2CK *restrict x_ck,
                      Mtrx3x3 *restrict g,
                      MtrxIdx3 a,
                      MtrxIdx3 b);

void prod_vuwdagger3x3(const Mtrx3x3 *restrict v,
                       const Mtrx3x3 *restrict u,
                       const Mtrx3x3 *restrict w,
                       Mtrx3x3 *restrict vuwdagger);

short power3x3Binomial(const Mtrx3x3 *restrict A,
                       const double omega,
                       Mtrx3x3 *restrict A_to_omega);

void updateSubLosAlamos(Mtrx3x3 *restrict w, Submtrx sub);

double inverse3x3(const Mtrx3x3 *restrict a, Mtrx3x3 *restrict a_inv);

short projectSU3(Mtrx3x3 *x);

void decomposeAlgebraSU3(const Mtrx3x3 *a, MtrxSU3Alg *a_components);

void projectSU3CabbiboMarinari(Mtrx3x3 *restrict w, Mtrx3x3 *restrict total_update);

void prod3x3_generic(Mtrx3x3 **restrict u, short number_of_matrices,
                     Mtrx3x3 *restrict result);

#endif  // SU3OPS_H