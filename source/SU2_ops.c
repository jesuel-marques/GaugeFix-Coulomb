#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>  //	Standard C header files

#include "math_ops.h"  //	Math operations
#include "lattice.h"
#include "SU2_ops.h"

// All matrices are in the Cayley-Klein representation
//	u=u[0] SU2_identity + i sum_i=1^3 u[i]sigma[i]
//	where sigma[i] are the Pauli matrices.

void print_matrix_2x2(const matrix_2x2 * restrict u, const char *name, const unsigned short decimal_places) {
    // Prints the matrix on screen with a given number of decimal places and 
    // adds a name on the top

    printf("\n\n %s \n", name);

    printf("{");
    for (SU2_color_idx  a = 0; a < 2; a++) {
        printf("{");

        for (SU2_color_idx  b = 0; b < 2; b++) {
            printf("%.*lf+I(%.*lf)", decimal_places, creal(u->m[ELM2x2(a, b)]), 
                                     decimal_places, cimag(u->m[ELM2x2(a, b)]));
            
            b != 2 -1 ?  printf(",") : 0 ;
        }

        a != 2 - 1 ?  printf("},\n") : 0 ;
    }

    printf("}}\n\n");

    getchar();
}

void copy_2x2(const mtrx_2x2_ck * restrict u, mtrx_2x2_ck * restrict u_copy) {
    // Copies u to u_copy

    for (SU2_color_idx a = 0; a < 4; a++) {

        u_copy -> m[a] = u -> m[a];

    }
}

inline void convert_from_ck(const mtrx_2x2_ck * restrict u_ck, matrix_2x2 * restrict u){
    
    u -> m[ELM2x2(0, 0)] =      u_ck -> m[0] 
                          + I * u_ck -> m[3];

    u -> m[ELM2x2(0, 1)] =      u_ck -> m[2]
                          + I * u_ck -> m[1];

    u -> m[ELM2x2(1, 0)] =   -  u_ck -> m[2] 
                          + I * u_ck -> m[1];

    u -> m[ELM2x2(1, 1)] =      u_ck -> m[0] 
                          - I * u_ck -> m[3];

}
// void set_null_2x2(mtrx_2x2_ck * restrict u) {
//     // Sets u to be the null Cayley-Klein 2x2 matrix

//     for (SU2_color_idx a = 0; a < 4; a++) {

//         u -> m[a] = 0.0;

//     }
// }

// void set_identity_2x2(mtrx_2x2_ck * restrict u) {
//     // Sets u to be the identity Cayley-Klein 2x2 matrix

//     u -> m[0] = 1.0;
//     u -> m[1] = 0.0;
//     u -> m[2] = 0.0;
//     u -> m[3] = 0.0;
// }

// void accumulate_2x2(const mtrx_2x2_ck * restrict u, mtrx_2x2_ck * restrict acc) {
//     // Accumulates the value of u into acc

//     for (SU2_color_idx a = 0; a < 4; a++) {

//         acc -> m[a] += u -> m[a];

//     }
// }

// void subtraction_2x2(const mtrx_2x2_ck * restrict u, 
//                      const mtrx_2x2_ck * restrict v, 
//                                 mtrx_2x2_ck * restrict u_minus_v) {
//     //  Calculates the difference between matrix u and matrix v
//     //  and returns result in u_minus_v

//     for (SU2_color_idx a = 0; a < 4; a++) {

//         u_minus_v -> m[a] = (u -> m[a]) 
//                           - (v -> m[a]);

//     }
// }

// work_data_type SU2_trace(const mtrx_2x2_ck * restrict u) {
//     //	Calculates trace of u

//     return 2.0 * (u -> m[0]);
//     //	In the Cayley-Klein representation, the trace is
//     //	twice the 0th component.
// }

inline work_data_type determinant_2x2(const mtrx_2x2_ck * restrict u) {
    //  Calculates the determinant of the matrix u

    work_data_type det_u = 0.0;

    for (SU2_color_idx i = 0; i <= 3; i++) {

        det_u += POW2(u -> m[i]);
        //	In the Cayley-Klein representation, the determinant
        //	is the sum of the squares of the components.

    }

    return det_u;
}

void herm_conj_2x2(const mtrx_2x2_ck * restrict u, mtrx_2x2_ck * restrict u_dagger) {
    // Calculates the hermitean conjugate to u
    // and returns result in u_dagger.

    u_dagger -> m[0] = u -> m[0];
    //	In the Cayley-Klein representation, the 0th
    //	component of the conjugate is the same...

    for (SU2_color_idx i = 1; i <= 3; i++) {

        u_dagger -> m[i] = -(u -> m[i]);
        //	And the 1, 2 and 3 components are the
        //	same up to a minus sign.

    }
}

void mult_scalar_2x2(const mtrx_2x2_ck * restrict u, const work_data_type alpha, 
                                                    mtrx_2x2_ck * restrict alpha_times_u) {
    //  Calculates multiplication of Cayley-Klein 2x2 matrix u by scalar alpha
    //  and returns result in alpha_times_u.

    for (SU2_color_idx a = 0; a < 4; a++) {

        alpha_times_u -> m[a] = alpha * (u -> m[a]);
        //	Mutiplying each entry.

    }
}

// static work_data_type SU2_inner_prod(const mtrx_2x2_ck * restrict u, const mtrx_2x2_ck * restrict v) {
//     //	Calculates the "scalar product" of two Cayley-Klein 2x2 matrices
//     //	in the Cayley-Klein representation.
//     //	Used in the product of two Cayley-Klein 2x2 matrices.

//     work_data_type inner_prod = (u -> m[0]) * (v -> m[0]);
//     //	The 0th component has a plus sign ...

//     for (SU2_color_idx b = 1; b < 4; b++) {

//         inner_prod += -(u -> m[b]) 
//                      * (v -> m[b]);

//     }
//     // and the 1, 2 and 3 components have a minus sign.

//     return inner_prod;
// }

// static void SU2_outer_product(const mtrx_2x2_ck * restrict u, 
//                               const mtrx_2x2_ck * restrict v, 
//                                         mtrx_2x2_ck * restrict outer_product) {
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

// void product_2x2(const mtrx_2x2_ck * restrict u, 
//                  const mtrx_2x2_ck * restrict v, 
//                             mtrx_2x2_ck * restrict uv) {
//     // Calculates product of 2 Cayley-Klein 2x2 matrices u e v
//     // and returns result in uv.

//     mtrx_2x2_ck u_cross_v;

//     uv -> m[0] = SU2_inner_prod(u, v);
//     //	In the Cayley-Klein representation, the 0th
//     //	components is given by the inner product...

//     SU2_outer_product(u, v, &u_cross_v);

//     for (SU2_color_idx a = 1; a <= 3; a++) {

//         uv -> m[a] = (u -> m[a]) * (v -> m[0]) 
//                    + (u -> m[0]) * (v -> m[a])
//                    - (u_cross_v.m[a]);
        
//     }
//     //	... and the 1, 2 e 3 components are given
//     //	by minus the cross product summed with
//     //	a term which mixes 0 and 1, 2, 3 components.
// }

// void product_three_2x2(const mtrx_2x2_ck * restrict u,
//                        const mtrx_2x2_ck * restrict v, 
//                        const mtrx_2x2_ck * restrict w, 
//                                 mtrx_2x2_ck * restrict uvw) {
//     //  Calculates product of 3 Cayley-Klein 2x2 matrices u, v and w
//     //  and returns result in uvw.

//     mtrx_2x2_ck uv;
  
//     product_2x2(u, v, &uv);
//     product_2x2(&uv, w, uvw);

// }

// void product_four_2x2(const mtrx_2x2_ck * restrict u, 
//                       const mtrx_2x2_ck * restrict v,
//                       const mtrx_2x2_ck * restrict w, 
//                       const mtrx_2x2_ck * restrict x, 
//                                 mtrx_2x2_ck * restrict uvwx) {
//     //  Calculates product of 4 Cayley-Klein 2x2 matrices u, v, w and x
//     //  and returns result in uvwx.

//     mtrx_2x2_ck uvw;

//     product_three_2x2(u, v, w, &uvw);
//     product_2x2(&uvw, x, uvwx);

// }

inline short SU2_projection(mtrx_2x2_ck * restrict u) {
    //	Projects matrix a to the group SU(2) returning SU(2) matrix a_SU2.

    mtrx_2x2_ck u_SU2;
    work_data_type det_u = determinant_2x2(u);
    if(det_u != 0){

        mult_scalar_2x2(u, 1.0 / sqrt(det_u), &u_SU2);
        copy_2x2(&u_SU2, u);
    
        return  0;

    }
    else{

        fprintf(stderr, "A matrix could not be projected to SU(2)");

        return -1;
    
    }
}