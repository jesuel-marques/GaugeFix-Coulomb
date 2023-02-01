/*
    SU3_ops performs 3x3 matrix operations.

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
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include <SU2_ops.h>
#include <SU3_ops.h>
#include <types.h>

#define DIGITS_MATRIX 18

/* Prints the 3x3 matrix to screen in Mathematica style. */
void printMatrix3x3(const Mtrx3x3 * restrict u) {
    /*
	 * Calls:
	 * =====
     * printf,
     * creal, cimag. 
     * 
     * Macros:
	 * ======
     * DIGITS_MATRIX,
     * Nc, ELEM_3X3, LOOP_3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:     3x3 matrix to be printed.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    printf("{");
    MtrxIdx3 a, b;
    
    LOOP_3(a) {
        printf("{");
        LOOP_3(b) {
            printf("%.*lf + I*(%.*lf)", DIGITS_MATRIX, creal(u -> m[ELEM_3X3(a, b)]), 
                                        DIGITS_MATRIX, cimag(u -> m[ELEM_3X3(a, b)]));
            
            b != Nc -1 ?  printf(", ") : 0 ;
        }
        a != Nc - 1 ?  printf("},\n") : 0 ;
    }
    printf("}}\n\n");
}

/* Copies a 3x3 matrix. */
void copy3x3(const Mtrx3x3 * restrict u, Mtrx3x3 * restrict u_copy) {
    /* 
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:         3x3 matrix to be copied,
     * Mtrx3x3 * u_copy:    copied matrix.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u_copy -> m[ELEM_3X3(a, b)] = u -> m[ELEM_3X3(a, b)];
    }
}

/* Sets a 3x3 matrix to be the 3x3 null matrix. */
inline void setNull3x3(Mtrx3x3 * restrict u) {
    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:         3x3 matrix to be set to null.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u -> m[ELEM_3X3(a, b)] = 0.0;
    }
}

/* Sets a 3x3 matrix to be the 3x3 identity matrix. */
inline void setIdentity3x3(Mtrx3x3 * restrict u) {
    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:         3x3 matrix to be set to the 3x3 identity matrix.
     * 
	 * Returns:
	 * =======
	 * 
	 */
    
    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u -> m[ELEM_3X3(a, b)] = a != b ? 0.0 : 1.0 ;
    }
}

/* Accumulates value of 3x3 matrix u in 3x3 matrix acc. */
inline void accumulate3x3(const Mtrx3x3 * restrict u, Mtrx3x3 * restrict acc) {
    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:         3x3 matrix to be accumulated,
     * Mtrx3x3 * acc:       accumulated 3x3 matrix.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        acc -> m[ELEM_3X3(a, b)] += u -> m[ELEM_3X3(a, b)];
    }
}

/* Subtracts two 3x3 matrices and put value in another one. */
void subtraction3x3(const Mtrx3x3 * restrict u, 
                     const Mtrx3x3 * restrict v,
                           Mtrx3x3 * restrict u_minus_v) {

    /* 
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:         3x3 matrix to be subtracted from,
     * Mtrx3x3 * v:         3x3 subtrahend matrix,
     * Mtrx3x3 * u_minus_v: difference matrix.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u_minus_v -> m[ELEM_3X3(a, b)] =   u -> m[ELEM_3X3(a, b)] 
                                         - v -> m[ELEM_3X3(a, b)];
    }
}

/* Obtains the trace of a 3x3 matrix. */
Scalar trace3x3(const Mtrx3x3 * restrict u) {

    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:         3x3 matrix whose trace is being calculated.
     * 
	 * Returns:
	 * =======
	 * Scalar containing the trace of the matrix, 
     * which is a complex number in double precision.
	 */

    return   u -> m[ELEM_3X3(0, 0)] 
           + u -> m[ELEM_3X3(1, 1)] 
           + u -> m[ELEM_3X3(2, 2)];
}

/* Obtains the determinant of a 3x3 matrix. */
inline Scalar determinant3x3(const Mtrx3x3 * restrict u) {

    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:     3x3 matrix whose determinant is being calculated.
     * 
	 * Returns:
	 * =======
	 * Scalar containing the determinant of the matrix, 
     * which is a complex number in double precision.
	 */

    //  The determinant is calculated in the usual manner  using Leibniz formula.
    return  u -> m[ELEM_3X3(0, 0)] * (  u->m[ELEM_3X3(1, 1)] * u->m[ELEM_3X3(2, 2)] 
                                      - u->m[ELEM_3X3(1, 2)] * u->m[ELEM_3X3(2, 1)]) 
          + u -> m[ELEM_3X3(0, 1)] * (  u->m[ELEM_3X3(1, 2)] * u->m[ELEM_3X3(2, 0)] 
                                      - u->m[ELEM_3X3(1, 0)] * u->m[ELEM_3X3(2, 2)]) 
          + u -> m[ELEM_3X3(0, 2)] * (  u->m[ELEM_3X3(1, 0)] * u->m[ELEM_3X3(2, 1)] 
                                      - u->m[ELEM_3X3(1, 1)] * u->m[ELEM_3X3(2, 0)]);
}

/* Calculates the hermitean conjugate of a 3x3 matrix. */
inline void hermConj3x3(const Mtrx3x3 * restrict u, 
                                Mtrx3x3 * restrict u_dagger) {

    /*
	 * Calls:
	 * =====
     * conj.
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:         3x3 matrix whose hermitean conjugate is being calculated,
     * Mtrx3x3 * u_dagger:  hermitean conjugated u.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    u_dagger -> m[ELEM_3X3(0, 0)] = conj(u -> m[ELEM_3X3(0, 0)]);
    u_dagger -> m[ELEM_3X3(1, 1)] = conj(u -> m[ELEM_3X3(1, 1)]);
    u_dagger -> m[ELEM_3X3(2, 2)] = conj(u -> m[ELEM_3X3(2, 2)]);

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u_dagger -> m[ELEM_3X3(b, a)] = conj(u -> m[ELEM_3X3(a, b)]);
        u_dagger -> m[ELEM_3X3(a, b)] = conj(u -> m[ELEM_3X3(b, a)]);
    }

}

/* Takes a 3x3 matrix u and calculates u - u^dagger - trace. */
inline void subtractionHermConjTrless3x3(const Mtrx3x3 * restrict u,
                                               Mtrx3x3 * restrict u_minus_udag_trless) {

    /*
	 * Calls:
	 * =====
     * creal, cimag, conj.
     * 
     * Macros:
	 * ======
     * Nc, ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:                    3x3 matrix whose traceless hermitean conjugated,
     *                                 subtraction is being calculated,
     * Mtrx3x3 * u_minus_udag_trless:  traceless hermitean conjugation subtraction u.
     * 
	 * Returns:
	 * =======
	 * 
	 */                              

    Scalar tr = 0.0;
    MtrxIdx3 a;

    LOOP_3(a){
        tr += u_minus_udag_trless -> m[ELEM_3X3(a, a)] 
            =     2 * I * cimag(u -> m[ELEM_3X3(a, a)]);
    }

    tr /= (double) Nc;

    LOOP_3(a){
        u_minus_udag_trless -> m[ELEM_3X3(a, a)] -= tr;
    }

    MtrxIdx3 b;
    LOOP_3X3(a, b) {
        if( a != b){
            u_minus_udag_trless -> m[ELEM_3X3(a, b)] =       u -> m[ELEM_3X3(a, b)] 
                                                      - conj(u -> m[ELEM_3X3(b, a)]);
            u_minus_udag_trless -> m[ELEM_3X3(b, a)] =       u -> m[ELEM_3X3(b, a)] 
                                                      - conj(u -> m[ELEM_3X3(a, b)]);
        }
    }

}

/* Calculates multiplication of 3x3 matrix u by Scalar alpha, which is a complex 
   number in double precision. */
inline void multByScalar3x3(const Scalar alpha, 
                            const Mtrx3x3 * restrict u,
                            Mtrx3x3 * restrict alpha_times_u) {

    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Scalar alpha:            scalar to multiply the matrix,
     * Mtrx3x3 * u:             3x3 matrix to be multiplied by alpha,
     * Mtrx3x3 * alpha_times_u: matrix resulting of multiplication of alpha by u.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        alpha_times_u -> m[ELEM_3X3(a, b)] = alpha * u -> m[ELEM_3X3(a, b)];
    }
}

/* Calculates multiplication of 3x3 matrix u by Scalar alpha, which is a complex 
   number in double precision. Result is returned in same matrix. */
inline void substMultScalar3x3(const Scalar alpha, 
                               Mtrx3x3 * restrict u) {

    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Scalar alpha:        scalar to multiply the matrix,
     * Mtrx3x3 * u:         3x3 matrix to be multiplied by alpha 
     *                      and in which the result will be put.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        u -> m[ELEM_3X3(a, b)] *= alpha;
    }
}

/* Calculates product of two 3x3 matrices. */
inline void prod3x3(const Mtrx3x3 * restrict u, 
                    const Mtrx3x3 * restrict v, 
                          Mtrx3x3 * restrict uv) {

    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:     first 3x3 matrix to be taken in the product,
     * Mtrx3x3 * v:     second 3x3 matrix to be taken in the product, 
     * Mtrx3x3 * uv:    result of product.
     *             
	 * Returns:
	 * =======
	 * 
	 */   
    
    MtrxIdx3 a, b, c;
    LOOP_3X3(a, b) {
        uv -> m[ELEM_3X3(a, b)] = 0.0;
        LOOP_3(c) {
            uv -> m[ELEM_3X3(a, b)] +=   (u -> m[ELEM_3X3(a, c)]) 
                                       * (v -> m[ELEM_3X3(c, b)]);
        }
    }
}

/* Calculates product of three 3x3 matrices. */
void prodThree3x3(const Mtrx3x3 * restrict u, 
                  const Mtrx3x3 * restrict v,
                  const Mtrx3x3 * restrict w,
                        Mtrx3x3 * restrict uvw) {

    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:     first 3x3 matrix to be taken in the product,
     * Mtrx3x3 * v:     second 3x3 matrix to be taken in the product,
     * Mtrx3x3 * w:     third 3x3 matrix to be taken in the product, 
     * Mtrx3x3 * uvw:   result of product.
     *             
	 * Returns:
	 * =======
	 * 
	 */

    MtrxIdx3 a, b, c, d;
    LOOP_3X3(a, b) {
        uvw -> m[ELEM_3X3(a, b)] = 0.0;
        LOOP_3(c) {
            LOOP_3(d) {
                uvw -> m[ELEM_3X3(a, b)] +=   (u -> m[ELEM_3X3(a, c)]) 
                                            * (v -> m[ELEM_3X3(c, d)]) 
                                            * (w -> m[ELEM_3X3(d, b)]);            
            }
        }
    }
}


/* Calculates product of four 3x3 matrices. */
void prodFour3x3(const Mtrx3x3 * restrict u, 
                 const Mtrx3x3 * restrict v,
                 const Mtrx3x3 * restrict w,
                 const Mtrx3x3 * restrict x, 
                       Mtrx3x3 * restrict uvwx) {

    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:     first 3x3 matrix to be taken in the product,
     * Mtrx3x3 * v:     second 3x3 matrix to be taken in the product,
     * Mtrx3x3 * w:     third 3x3 matrix to be taken in the product,
     * Mtrx3x3 * x:     fourth 3x3 matrix to be taken in the product,
     * Mtrx3x3 * uvwx:  result of product.
     *             
	 * Returns:
	 * =======
	 * 
	 */

    MtrxIdx3 a, b, c, d, e;
    LOOP_3X3(a, b) {
        uvwx -> m[ELEM_3X3(a, b)] = 0.0;
        LOOP_3(c) {
            LOOP_3(d) {
                LOOP_3(e) {
                    uvwx -> m[ELEM_3X3(a, b)] +=   (u -> m[ELEM_3X3(a, c)]) 
                                                 * (v -> m[ELEM_3X3(c, d)]) 
                                                 * (w -> m[ELEM_3X3(d, e)]) 
                                                 * (x -> m[ELEM_3X3(e, b)]);
                }
            }
        }
    }
}

/* Calculates matrix product between g and acc_prod and 
   accumulates result in acc_prod. */
inline void accumLeftProd3x3(const Mtrx3x3 * restrict g,
                                   Mtrx3x3 * restrict acc_prod) {

    /* 
	 * Calls:
	 * =====
     * copy3x3, prod3x3.
     * 
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * g:         3x3 matrix to be left product accumulated in the acc_prod
     * Mtrx3x3 * acc_prod:  3x3 matrix which holds the accumulated product
     *             
	 * Returns:
	 * =======
	 * 
	 */

    Mtrx3x3 aux_prod;
 
    copy3x3(acc_prod, &aux_prod);
    prod3x3(g, &aux_prod, acc_prod);

}

/* Calculates matrix product between acc_prod and g and accumulates result in acc_prod. 
   */
inline void accumRightProd3x3(      Mtrx3x3 * restrict acc_prod, 
                              const Mtrx3x3 * restrict g) {

    /*
	 * Calls:
	 * =====
     * copy3x3, prod3x3.
     * 
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
     * Mtrx3x3 * acc_prod:  3x3 matrix which holds the accumulated product
	 * Mtrx3x3 * g:         3x3 matrix to be right product accumulated in the acc_prod
     *             
	 * Returns:
	 * =======
	 * 
	 */
    
    Mtrx3x3 aux_prod;

    copy3x3(acc_prod, &aux_prod);
    prod3x3(&aux_prod, g, acc_prod);

}

/* Accumulates in g the product of a SU(2) Cabbibo-Marinari submatrix in
   Cayley-Klein form x_ck with g.  */
inline void accumProdSU2_3x3(const Mtrx2x2CK * restrict x_ck, 
                                   Mtrx3x3 * restrict g, 
                                   MtrxIdx3 a, 
                                   MtrxIdx3 b ) {
    /* 
	 * Calls:
	 * =====
     * convertFromCK.
     * 
     * Macros:
	 * ======
     * ELM_2X2,
     * LOOP_3, ELEM_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
     * Mtrx2x2CK * restrict x_ck:   SU(2) matrix to be multiplied with 3x3 matrix,
     * Mtrx3x3 * restrict g:        3x3 matrix to be multplied from left with SU(2)
     *                              with,
     * MtrxIdx3 a, b:               the relevant line and column of the 3x3 matrix.
     *             
	 * Returns:
	 * =======
	 * 
	 */

    Scalar xg1, xg2;
    //  Auxiliary variables
    
    Mtrx2x2 x;
    //  the Cayley-Klein matrix will be first converted to a regular 2x2 matrix
    convertFromCK(x_ck, &x);

    MtrxIdx3 c;
    LOOP_3(c) {
        //  Modifying only the terms in g that will actually receive contributions

        xg1 =   x.m[ELM_2X2(0, 0)] * g -> m[ELEM_3X3(a, c)] 
              + x.m[ELM_2X2(0, 1)] * g -> m[ELEM_3X3(b, c)];
        xg2 =   x.m[ELM_2X2(1, 0)] * g -> m[ELEM_3X3(a, c)] 
              + x.m[ELM_2X2(1, 1)] * g -> m[ELEM_3X3(b, c)];

        g -> m[ELEM_3X3(a, c)] = xg1;
        g -> m[ELEM_3X3(b, c)] = xg2;
    }

}

/* Calculates matrix product between v, u and w dagger. */
inline void prod_vuwdagger3x3(const Mtrx3x3 * restrict v, 
                              const Mtrx3x3 * restrict u, 
                              const Mtrx3x3 * restrict w, 
                                    Mtrx3x3 * restrict vuwdagger ) {

    /*
	 * Calls:
	 * =====
     * conj
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * u:             first 3x3 matrix to be taken in the product,
     * Mtrx3x3 * v:             middle 3x3 matrix to be taken in the product,
     * Mtrx3x3 * w:             3x3 matrix to be taken the hermitean conjugate and then, 
     *                          taken in the product,
     * Mtrx3x3 * vuwdagger:     result of product.
     *             
	 * Returns:
	 * =======
	 * 
	 */
    
    MtrxIdx3 a, b, c, d;
    LOOP_3X3(a, b) {
        vuwdagger -> m[ELEM_3X3(a, b)] = 0.0;
        LOOP_3(c) {
            LOOP_3(d) {
                vuwdagger -> m[ELEM_3X3(a, b)] +=       (v -> m[ELEM_3X3(a, c)]) 
                                                  *     (u -> m[ELEM_3X3(c, d)]) 
                                                  * conj(w -> m[ELEM_3X3(b, d)]);
                                                //  wdagger_db is conj(w_bd)
            }
        }
    }
}

/* Calculates A^omega according to the formula (which works for A close to identity)
   A^omega = Proj_SU3(I + omega * (A-I) + ...). */
inline short power3x3Binomial(const Mtrx3x3 * restrict A, 
                              const double omega, 
                              Mtrx3x3 * restrict A_to_omega) {
    /*
	 * Calls:
	 * =====
     * projectSU3
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * A:             3x3 matrix whose power is being taken
     * double omega:            the power that is being taken
     * Mtrx3x3 * A_to_omega:    result of the exponentiation
     *             
	 * Returns:
	 * =======
	 * 
	 */

    //  A^omega = Proj_SU3(I + omega * (A-I) + ...)
    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {

        A_to_omega -> m[ELEM_3X3(a, b)] = a == b ? 
                                          1.0 + omega * (A -> m[ELEM_3X3(a, b)] - 1.0):
                                                omega * (A -> m[ELEM_3X3(a, b)]      );
            
    }

    return projectSU3(A_to_omega);
}

/* Performs a Los Alamos update to w with submatrix indicated by sub. The Los Alamos 
   update is a SU(2) Cabbibo-Marinari submatrix of SU(3) which maximizes the value of
   ReTr(w). */
inline void updateSubLosAlamos(Mtrx3x3 * restrict w, const Submtrx sub) {

    /*  
	 * Calls:
	 * =====
     * creal, cimag,
     * projectSU2,
     * accumProdSU2_3x3.
     * 
	 * Macros:
	 * ======
     * ELEM_3X3.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * w:     matrix to be updated; comes from accumulating links,
     * Submtrx sub:     index of SU(2) Cabbibo-Marinari submatrix of SU(3).
     *  
	 * Returns:
	 * =======
	 *
	 */

    MtrxIdx3 a, b;

    a = sub == T ? 1 : 0;
    b = sub == R ? 1 : 2;

    /* a and b will be the line and colums index
    for each Cabbibo-Marinari matrix,  */

    Mtrx2x2CK LA_sub_update;

    LA_sub_update.m[0] =  creal( w -> m[ELEM_3X3(a, a)]
                               + w -> m[ELEM_3X3(b, b)]);
    LA_sub_update.m[1] = -cimag( w -> m[ELEM_3X3(a, b)] 
                               + w -> m[ELEM_3X3(b, a)]);
    LA_sub_update.m[2] = -creal( w -> m[ELEM_3X3(a, b)]
                               - w -> m[ELEM_3X3(b, a)]);
    LA_sub_update.m[3] = -cimag( w -> m[ELEM_3X3(a, a)] 
                               - w -> m[ELEM_3X3(b, b)]);
  

    /*  The SU(2) matrix corresponding to the Cabbibo-Marinari submatrix is built
        from w according to the formulae above. This is what maximizes the functional 
        locally for each submatrix. This can be proven using maximization with 
        constraints, where the constraint is that mtrx_SU2 should have determinant equal
        to one, when written in the Caylel-Klein form, which guarantees that it is
         unitary. */

    if(!projectSU2(&LA_sub_update)) {

        /* If projectSU2 is succesful, it will return 0. Then w can be updated */
        accumProdSU2_3x3(&LA_sub_update, w, a, b);
        return;
    }
    else{
        /*  If projectSU2 is unsuccesful, this means that mtrx_SU2 was 0 and w remains
            what it was. The program will then carry on to the next submatrix update. */
        return;
    }
}

/* Calculates the inverse of a 3 by 3 matrix x explicitly. */
inline double inverse3x3(const Mtrx3x3 * restrict a,
                               Mtrx3x3 * restrict a_inv) {
                                    
    /*
	 * Calls:
	 * =====
     * determinant3x3.
     * 
     * Macros:
	 * ======
     * ELEM_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * a:         3x3 matrix whose inverse is wanted,
     * Mtrx3x3 * a_inv:     inverse of a if exists, if it doesn't, the null matrix.
     * 
	 * Returns:
	 * =======
	 * The determinant of a in double precision, can be used to check if the determinant
     * of matrix was calculated.
	 */
	
	Scalar a_det = determinant3x3(a);
    
    if(a_det != 0.0) {

        a_inv->m[ELEM_3X3(0, 0)] =  (a->m[ELEM_3X3(1, 1)] * a->m[ELEM_3X3(2, 2)] 
                                   - a->m[ELEM_3X3(1, 2)] * a->m[ELEM_3X3(2, 1)])/a_det;
        a_inv->m[ELEM_3X3(0, 1)] =  (a->m[ELEM_3X3(0, 2)] * a->m[ELEM_3X3(2, 1)] 
                                   - a->m[ELEM_3X3(0, 1)] * a->m[ELEM_3X3(2, 2)])/a_det;
        a_inv->m[ELEM_3X3(0, 2)] =  (a->m[ELEM_3X3(0, 1)] * a->m[ELEM_3X3(1, 2)] 
                                   - a->m[ELEM_3X3(0, 2)] * a->m[ELEM_3X3(1, 1)])/a_det;

        a_inv->m[ELEM_3X3(1, 0)] =  (a->m[ELEM_3X3(1, 2)] * a->m[ELEM_3X3(2, 0)] 
                                   - a->m[ELEM_3X3(1, 0)] * a->m[ELEM_3X3(2, 2)])/a_det;
        a_inv->m[ELEM_3X3(1, 1)] =  (a->m[ELEM_3X3(0, 0)] * a->m[ELEM_3X3(2, 2)] 
                                   - a->m[ELEM_3X3(0, 2)] * a->m[ELEM_3X3(2, 0)])/a_det;
        a_inv->m[ELEM_3X3(1, 2)] =  (a->m[ELEM_3X3(0, 2)] * a->m[ELEM_3X3(1, 0)] 
                                   - a->m[ELEM_3X3(0, 0)] * a->m[ELEM_3X3(1, 2)])/a_det;

        a_inv->m[ELEM_3X3(2, 0)] =  (a->m[ELEM_3X3(1, 0)] * a->m[ELEM_3X3(2, 1)] 
                                   - a->m[ELEM_3X3(1, 1)] * a->m[ELEM_3X3(2, 0)])/a_det;
        a_inv->m[ELEM_3X3(2, 1)] =  (a->m[ELEM_3X3(0, 1)] * a->m[ELEM_3X3(2, 0)] 
                                   - a->m[ELEM_3X3(0, 0)] * a->m[ELEM_3X3(2, 1)])/a_det;
        a_inv->m[ELEM_3X3(2, 2)] =  (a->m[ELEM_3X3(0, 0)] * a->m[ELEM_3X3(1, 1)] 
                                   - a->m[ELEM_3X3(0, 1)] * a->m[ELEM_3X3(1, 0)])/a_det;
    
    }
    else{
        setNull3x3(a_inv);
    }
    
    return a_det;
}

/* Takes a color vector and normalizes it. */
inline static int normalize3Vec(Vec3 * restrict v) {
    /*
	 * Calls:
	 * =====
     * fabs, sqrt.
     * 
     * Macros:
	 * ======
     * POW2, 
     * LOOP_3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * v:     color vector to be normalized.
     * 
	 * Returns:
	 * =======
	 * 0 if successful, -1 otherwise. Normalization will not be successful if vector was
     * the null vector.
	 */

    double r = 0.0;

    MtrxIdx3 a;

    LOOP_3(a) {
        //  Calculates the sum of the absolute
        //  value squared of each component
        r += (fabs(v -> m[a])) * (fabs(v -> m[a]));

    }
    if(r != 0) {
        //  If vector is not null, then normalize it
        r = 1.0 / sqrt(r);

        LOOP_3(a) {

            v -> m[a] *= r;

        }
        return 0;
    }
    else{
        return -1;
    }
}

/* Calculates the conjugate of the cross product of color vectors u and v. */
inline static void crossproductConj(Vec3 u, 
                                    Vec3 v,
                                    Vec3* u_cross_v_star) {

    /*
	 * Calls:
	 * =====
     * conj.
     * 
     * Macros:
	 * ======
     * LOOP_3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
     * Mtrx3x3 * u:     the first vector in the cross product,
	 * Mtrx3x3 * v:     the second vector in the cross product.
     * 
	 * Returns:
	 * =======
	 * 
     */

    u_cross_v_star -> m[0] = conj(u.m[1] * v.m[2] 
                                - u.m[2] * v.m[1]);

    u_cross_v_star -> m[1] = conj(u.m[2] * v.m[0] 
                                - u.m[0] * v.m[2]);

    u_cross_v_star -> m[2] = conj(u.m[0] * v.m[1] 
                                - u.m[1] * v.m[0]);

}

/* Projects matrix x to the group SU(3) returning SU(3) matrix in x at the end.
   Follows method found in openqcd fastsum code.
   https://gitlab.com/fastsum/openqcd-fastsum/-/tree/master/ */
inline short projectSU3(Mtrx3x3 * restrict x) {

    /*
	 * Calls:
	 * =====
     * normalize3Vec, crossproductConj.
     * 
     * Macros:
	 * ======
     * ELEM_3X3, LOOP_3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
     * Mtrx3x3 * x:     the 3x3 matrix to be projected to SU(3).
     * 
	 * Returns:
	 * =======
	 * 
     */

    short exit_status;

    Vec3 v0, v1, v2;

    MtrxIdx3 b;
    LOOP_3(b) {

        v0.m[b] = x -> m[ELEM_3X3(0, b)];
        v1.m[b] = x -> m[ELEM_3X3(1, b)];

    }

    exit_status = normalize3Vec(&v0);
    crossproductConj(v0, v1, &v2);

    exit_status |= normalize3Vec(&v2);
    crossproductConj(v2, v0, &v1);

    LOOP_3(b) {

        x -> m[ELEM_3X3(0, b)] = v0.m[b];
        x -> m[ELEM_3X3(1, b)] = v1.m[b];
        x -> m[ELEM_3X3(2, b)] = v2.m[b];

    }
        
    if(exit_status == 0) {
        return 0;
    }
    else{
        return -1;
    }
}

/* Decomposes A in the components of the alfebra of SU(3) using the 
   Gell-Mann matrices as a basis and following the conventions of Gattringer 
   (https://link.springer.com/book/10.1007/978-3-642-01850-3). 
   Formulas for the coefficients obtained in Mathematica. */
void decomposeAlgebraSU3(const Mtrx3x3    * restrict a, 
                               MtrxSU3Alg * restrict a_components) {
    
    /*
	 * Calls:
	 * =====
     * sqrt,
     * creal, cimag.
     * 
     * Macros:
	 * ======
     * ELEM_3X3.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
     * Mtrx3x3 * a:             the SU(3) algebra matrix whose components are wanted,
	 * Mtrx3x3 * a_components:  the Gell-Mann basis components of the a matrix.
     * 
	 * Returns:
	 * =======
	 * 
     */

    a_components -> m[0] =       0.0;

    a_components -> m[1] =       ( creal(a -> m[ELEM_3X3(0, 1)]) 
                                 + creal(a -> m[ELEM_3X3(1, 0)]));

    a_components -> m[2] =       (-cimag(a -> m[ELEM_3X3(0, 1)]) 
                                 + cimag(a -> m[ELEM_3X3(1, 0)]));

    a_components -> m[3] =       ( creal(a -> m[ELEM_3X3(0, 0)]) 
                                 - creal(a -> m[ELEM_3X3(1, 1)]));

    a_components -> m[4] =       ( creal(a -> m[ELEM_3X3(0, 2)]) 
                                 + creal(a -> m[ELEM_3X3(2, 0)]));

    a_components -> m[5] =       (-cimag(a -> m[ELEM_3X3(0, 2)]) 
                                 + cimag(a -> m[ELEM_3X3(2, 0)]));

    a_components -> m[6] =       ( creal(a -> m[ELEM_3X3(1, 2)]) 
                                 + creal(a -> m[ELEM_3X3(2, 1)]));

    a_components -> m[7] =       (-cimag(a -> m[ELEM_3X3(1, 2)]) 
                                 + cimag(a -> m[ELEM_3X3(2, 1)]));

    a_components -> m[8] =       ( creal(a -> m[ELEM_3X3(0, 0)]) 
                                 + creal(a -> m[ELEM_3X3(1, 1)]) 
                           - 2.0 * creal(a -> m[ELEM_3X3(2, 2)])) / sqrt(3.0);
}

#define CABBIBO_MARINARI_HITS 300

/* Alternative projection method, following */
void projectSU3CabbiboMarinari(Mtrx3x3 * restrict w, 
                               Mtrx3x3 * restrict w_projected) {
    
    /*  
	 * Calls:
	 * =====
     * inverse3x3, updateSubLosAlamos, setIdentity3x3, prod3x3, hermConj3x3.
     * 
     * Macros:
	 * ======
     * CABBIBO_MARINARI_HITS.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
     * Mtrx3x3 * w:             the SU(3) algebra matrix whose components are wanted,
	 * Mtrx3x3 * a_components:  the Gell-Mann basis components of the a matrix.
     * 
	 * Returns:
	 * =======
	 * 
     */

    Mtrx3x3 w_inv_old;
    Mtrx3x3 total_update_conj;
    /* Calculates the inverse of w in the beginning. The program will update w 
       successively and to extract what was the combined update, we can multiply from
       the right by the old inverse. */

    if(inverse3x3(w, &w_inv_old)) {

        /* Local maximization is attained iteratively in SU(3),
           thus we need to make many hits ... */
        for(unsigned short hits = 1; hits <= CABBIBO_MARINARI_HITS; hits++) {

            /* ... and each hit contains the Cabbibo-Marinari subdivision */
            for(Submtrx sub = R; sub <= T; sub++) {
                /* Submatrices are indicated by numbers from 0 to 2 with 
                   codenames R, S and T */

                updateSubLosAlamos(w, sub);
            }
        }
    }
    else{
        //  if w has no inverse, update will be given by the identity
        setIdentity3x3(w_projected);
        return;
    }
    
    prod3x3(w, &w_inv_old, &total_update_conj);
    hermConj3x3(&total_update_conj, w_projected);

    //	Updates matrix to total_update. It is the accumulated updates from the hits.
}