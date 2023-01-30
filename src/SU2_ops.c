/*
    SU2_ops performs 2x2 matrix operations.

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
#include <tgmath.h>

#include <misc.h>
#include <SU2_ops.h>
#include <types.h>


// All matrices are in the Cayley-Klein representation
//	u=u[0] SU2_identity + i sum_i=1^3 u[i]sigma[i]
//	where sigma[i] are the Pauli matrices.

// void printMtrx2x2(const Mtrx2x2 * restrict u, 
//                     const char *name, 
//                     const unsigned short decimal_places) {
//     // Prints the matrix on screen with a given number of decimal places and 
//     // adds a name on the top

//     printf("\n\n %s \n", name);

//     printf("{");
//     for(MtrxIdx2  a = 0; a < 2; a++) {
//         printf("{");

//         for(MtrxIdx2  b = 0; b < 2; b++) {
//             printf("%.*lf + I(%.*lf)", decimal_places, creal(u -> m[ELM_2X2(a, b)]), 
//                                        decimal_places, cimag(u -> m[ELM_2X2(a, b)]));
            
//             b != 2 -1 ?  printf(",") : 0 ;
//         }

//         a != 2 - 1 ?  printf("},\n") : 0 ;
//     }

//     printf("}}\n\n");

//     getchar();
// }

/* Copies a 2x2 matrix in the Cayley-Klein representation. */
void copy2x2(const Mtrx2x2CK * restrict u, 
             Mtrx2x2CK * restrict u_copy) {

    /* 
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * LOOP_2_CK.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx2x2CK * u:         2x2 Cayley-Klein matrix to be copied,
     * Mtrx2x2CK * u_copy:    copied matrix.
     * 
	 * Returns:
	 * =======
	 * 
	 */
    
    MtrxIdx2 a;
    LOOP_2_CK(a) {

        u_copy -> m[a] = u -> m[a];

    }
}

/* Converts a 2x2 matrix from the Cayley-Klein representation to a 2x2 usual matrix. */
inline void convertFromCK(const Mtrx2x2CK * restrict u_ck, 
                          Mtrx2x2 * restrict u) {

    /* 
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * ELM_2X2.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx2x2CK * u_ck:    2x2 Cayley-Klein matrix to be converted,
     * Mtrx2x2 * u:         2x2 usual matrix.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    u -> m[ELM_2X2(0, 0)] =      u_ck -> m[0] 
                           + I * u_ck -> m[3];

    u -> m[ELM_2X2(0, 1)] =      u_ck -> m[2]
                           + I * u_ck -> m[1];

    u -> m[ELM_2X2(1, 0)] =   -  u_ck -> m[2] 
                           + I * u_ck -> m[1];

    u -> m[ELM_2X2(1, 1)] =      u_ck -> m[0] 
                           - I * u_ck -> m[3];

}


// void setNull2x2(Mtrx2x2CK * restrict u) {
//     // Sets u to be the null Cayley-Klein 2x2 matrix
//     for(MtrxIdx2 a = 0; a < 4; a++) {

//         u -> m[a] = 0.0;

//     }
// }


// void setIdentity2x2(Mtrx2x2CK * restrict u) {
//     // Sets u to be the identity Cayley-Klein 2x2 matrix
//     u -> m[0] = 1.0;
//     u -> m[1] = 0.0;
//     u -> m[2] = 0.0;
//     u -> m[3] = 0.0;
// }


// void accumulate2x2(const Mtrx2x2CK * restrict u, 
//                           Mtrx2x2CK * restrict acc) {
//     // Accumulates the value of u into acc
//     for(MtrxIdx2 a = 0; a < 4; a++) {

//         acc -> m[a] += u -> m[a];

//     }
// }


// void subtraction2x2(const Mtrx2x2CK * restrict u, 
//                      const Mtrx2x2CK * restrict v, 
//                            Mtrx2x2CK * restrict u_minus_v) {
//     //  Calculates the difference between matrix u and matrix v
//     //  and returns result in u_minus_v

//     for(MtrxIdx2 a = 0; a < 4; a++) {

//         u_minus_v -> m[a] = (u -> m[a]) 
//                           - (v -> m[a]);

//     }
// }


// Scalar SU2Trace(const Mtrx2x2CK * restrict u) {
//     //	Calculates trace of u
//     return 2.0 * (u -> m[0]);
//     //	In the Cayley-Klein representation, the trace is
//     //	twice the 0th component.
// }

/* Obtains the determinant of a 2x2 Cayley-Klein matrix. */
inline Scalar determinant2x2(const Mtrx2x2CK * restrict u) {
    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * POW2,
     * LOOP_2_CK.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Mtrx2x2CK * u:     2x2 Cayley-Klein matrix whose determinant is being calculated.
     * 
	 * Returns:
	 * =======
	 * Scalar containing the determinant of the matrix, which is a real number cast 
     * to a complex number in double precision.
	 */
    
    Scalar det_u = 0.0;
    MtrxIdx2 a;
    LOOP_2_CK(a) {

        det_u += POW2(u -> m[a]);
        //	In the Cayley-Klein representation, the determinant
        //	is the sum of the squares of the components.

    }

    return det_u;
}


// void hermConj2x2(const Mtrx2x2CK * restrict u, 
//                          Mtrx2x2CK * restrict u_dagger) {
//     // Calculates the hermitean conjugate to u
//     // and returns result in u_dagger.
//     u_dagger -> m[0] = u -> m[0];
//     //	In the Cayley-Klein representation, the 0th
//     //	component of the conjugate is the same...
//     MtrxIdx2 a;
//     LOOP_2_CK_i(a) {

//         u_dagger -> m[a] = -(u -> m[a]);
//         //	And the 1, 2 and 3 components are the
//         //	same up to a minus sign.

//     }
// }

/* Calculates multiplication of 2x2 Cayley-Klein matrix u by Scalar alpha, 
   which is a complex number in double precision. */
void multByScalar2x2(const Scalar alpha,
                     const Mtrx2x2CK * restrict u, 
                     Mtrx2x2CK * restrict alpha_times_u) {
    
    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * LOOP_2_CK.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
	 * Scalar alpha:                scalar to multiply the matrix,
     * Mtrx2x2CK * u:               2x2 Cayley-Klein matrix to be multiplied by alpha,
     * Mtrx2x2CK * alpha_times_u:   matrix resulting of multiplication of alpha by u.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    MtrxIdx2 a;
    LOOP_2_CK(a) {
        alpha_times_u -> m[a] = alpha * (u -> m[a]);
    }
}


// static Scalar SU2InnerProd(const Mtrx2x2CK * restrict u, 
//                              const Mtrx2x2CK * restrict v) {
//     //	Calculates the "Scalar product" of two Cayley-Klein 2x2 matrices
//     //	in the Cayley-Klein representation.
//     //	Used in the product of two Cayley-Klein 2x2 matrices.
//     Scalar inner_prod = (u -> m[0]) * (v -> m[0]);
//     //	The 0th component has a plus sign ...

//     for(MtrxIdx2 b = 1; b < 4; b++) {

//         inner_prod += -(u -> m[b]) 
//                      * (v -> m[b]);

//     }
//     // and the 1, 2 and 3 components have a minus sign.

//     return inner_prod;
// }


// static void SU2OuterProduct(const Mtrx2x2CK * restrict u, 
//                               const Mtrx2x2CK * restrict v, 
//                                     Mtrx2x2CK * restrict outer_product) {
//     //	Calculates the "outer product" of two Cayley-Klein 2x2 matrices
//     //	in the Cayley-Klein representation, and returns result in outer_product
//     //	Used in the product of two Cayley-Klein 2x2 matrices.

//     //	Actually, product taken from the 1, 2 and 3 components only,
//     //	in the usual way, as in regular linear algebra.
//     outer_product -> m[0] = 0.0;
//     outer_product -> m[1] = (u -> m[2]) * (v -> m[3]) 
//                           - (u -> m[3]) * (v -> m[2]);
//     outer_product -> m[2] = (u -> m[3]) * (v -> m[1]) 
//                           - (u -> m[1]) * (v -> m[3]);
//     outer_product -> m[3] = (u -> m[1]) * (v -> m[2]) 
//                           - (u -> m[2]) * (v -> m[1]);
// }


// void product2x2(const Mtrx2x2CK * restrict u, 
//                  const Mtrx2x2CK * restrict v, 
//                        Mtrx2x2CK * restrict uv) {
//     // Calculates product of 2 Cayley-Klein 2x2 matrices u e v
//     // and returns result in uv.
//     Mtrx2x2CK u_cross_v;

//     uv -> m[0] = SU2InnerProd(u, v);
//     //	In the Cayley-Klein representation, the 0th
//     //	components is given by the inner product...

//     SU2OuterProduct(u, v, &u_cross_v);

//     for(MtrxIdx2 a = 1; a <= 3; a++) {

//         uv -> m[a] = (u -> m[a]) * (v -> m[0]) 
//                    + (u -> m[0]) * (v -> m[a])
//                    - (u_cross_v.m[a]);
        
//     }
//     //	... and the 1, 2 e 3 components are given
//     //	by minus the cross product summed with
//     //	a term which mixes 0 and 1, 2, 3 components.
// }


// void productThree2x2(const Mtrx2x2CK * restrict u,
//                        const Mtrx2x2CK * restrict v, 
//                        const Mtrx2x2CK * restrict w, 
//                              Mtrx2x2CK * restrict uvw) {
//     //  Calculates product of 3 Cayley-Klein 2x2 matrices u, v and w
//     //  and returns result in uvw.
//     Mtrx2x2CK uv;
  
//     product2x2(u, v, &uv);
//     product2x2(&uv, w, uvw);

// }


// void productFour2x2(const Mtrx2x2CK * restrict u, 
//                       const Mtrx2x2CK * restrict v,
//                       const Mtrx2x2CK * restrict w, 
//                       const Mtrx2x2CK * restrict x, 
//                             Mtrx2x2CK * restrict uvwx) {
//     //  Calculates product of 4 Cayley-Klein 2x2 matrices u, v, w and x
//     //  and returns result in uvwx.
//     Mtrx2x2CK uvw;

//     productThree2x2(u, v, w, &uvw);
//     product2x2(&uvw, x, uvwx);

// }

/* Projects 2x2 Cayley-Klein matrix u to the group SU(2) returning SU(2) matrix in
   u at the end. */
inline short projectSU2(Mtrx2x2CK * restrict u) {

    /*
	 * Calls:
	 * =====
     * copy2x2, determinant2x2, multByScalar2x2.
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
     * Mtrx3x3 * u:     the 2x2 Cayley-Klein matrix to be projected to SU(2).
     * 
	 * Returns:
	 * =======
	 * 
     */

    Mtrx2x2CK u_SU2;
    Scalar det_u = determinant2x2(u);
    if(det_u != 0) {
        multByScalar2x2(1.0 / sqrt(det_u), u, &u_SU2);
        copy2x2(&u_SU2, u);
        return  0;
    }
    else{
        return -1;
    }
}