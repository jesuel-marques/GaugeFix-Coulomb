/*
    Header to SU2_ops.c, which performs 2x2 matrix operations.

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

#ifndef SU2OPS_H
#define SU2OPS_H

#include <types.h>

#define M_PI 3.14159265358979323846 /* pi */

// 2x2 types

/* 2x2 matrix index */
typedef unsigned short MtrxIdx2;

/* 2x2 matrices */
typedef struct {
    Scalar m[2 * 2];
} Mtrx2x2;

/*  2x2 matrices in the Cayley-Klein form
    m = m[0] SU2_identity + i sum_i=1^3 m[i]sigma[i] */
typedef struct {
    double m[4];
} Mtrx2x2CK;

/* used to get the elements of 2x2 matrices */
#define ELM_2X2(a, b) (a) * 2 + (b)

/* used to loop through 2x2 indices for matrices in Cayley-Klein form */
#define LOOP_2_CK(a) for (a = 0; a < 4; a++)

/* used to loop through 2x2 indices for matrices in Cayley-Klein form
   excluding the 0th element */
#define LOOP_2_CK_i(a) for (a = 1; a < 4; a++)

/* used to loop through color indices for 2x2 matrices */
#define LOOP_2X2(a, b)      \
    for (a = 0; a < 2; a++) \
        for (b = 0; b < 2; b++)

//  void printMtrx2x2(const Mtrx2x2 * restrict u,
//                     const char *name,
//                     const unsigned short decimal_places);

void copy2x2(const Mtrx2x2CK* restrict u, Mtrx2x2CK* restrict u_copy);

void convertFromCK(const Mtrx2x2CK* restrict u_ck, Mtrx2x2* restrict u);

// void setNull2x2(Mtrx2x2CK * u), setIdentity2x2(Mtrx2x2CK * u);

// void accumulate2x2(const Mtrx2x2CK * u, Mtrx2x2CK* acc);
// void subtraction2x2(const Mtrx2x2CK* u, const Mtrx2x2CK* v, Mtrx2x2CK* u_minus_v);

// Scalar SU2Trace(const Mtrx2x2CK *u);

Scalar determinant2x2(const Mtrx2x2CK* u);
void hermConj2x2(const Mtrx2x2CK* u, Mtrx2x2CK* u_dagger);

void multByScalar2x2(const Scalar alpha,
                     const Mtrx2x2CK* restrict u,
                     Mtrx2x2CK* restrict alpha_times_u);

void product2x2(const Mtrx2x2CK* u, const Mtrx2x2CK* v, Mtrx2x2CK* uv);
// productThree2x2(const Mtrx2x2CK* u,
//                   const Mtrx2x2CK* v,
//                   const Mtrx2x2CK* w,
//                         Mtrx2x2CK* uvw),
// productFour2x2(const Mtrx2x2CK* u,
//                  const Mtrx2x2CK* v,
//                  const Mtrx2x2CK* w,
//                  const Mtrx2x2CK* x,
//                        Mtrx2x2CK* uvwx);

Scalar projectSU2(Mtrx2x2CK* restrict u);

#endif  // SU2OPS_H