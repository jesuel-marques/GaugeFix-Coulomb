#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>  //	Standard C header files

#include "math_ops.h"  //	Math operations
#include "lattice.h"
#include "SU2_ops.h"

// All matrices are in the Cayley-Klein representation
//	u=u[0] SU2_identity + i sum_i=1^3 u[i]sigma[i]
//	where sigma[i] are the Pauli matrices.

void copy_2x2(const matrix_2x2_ck *u, matrix_2x2_ck *u_copy) {
    // Copies u to u_copy

    for (SU2_color_index a = 0; a < 4; a++) {

        u_copy -> m[a] = u -> m[a];

    }
}

void set_to_null_2x2(matrix_2x2_ck *u) {
    // Sets u to be the null Cayley-Klein 2x2 matrix

    for (SU2_color_index a = 0; a < 4; a++) {

        u -> m[a] = 0.0;

    }
}

void set_to_identity_2x2(matrix_2x2_ck *u) {
    // Sets u to be the identity Cayley-Klein 2x2 matrix

    u -> m[0] = 1.0;
    u -> m[1] = 0.0;
    u -> m[2] = 0.0;
    u -> m[3] = 0.0;
}

void accumulate_2x2(const matrix_2x2_ck *u, matrix_2x2_ck *acc) {
    // Accumulates the value of u into acc

    for (SU2_color_index a = 0; a < 4; a++) {

        acc -> m[a] += u -> m[a];

    }
}

void subtraction_2x2(const matrix_2x2_ck *u, const matrix_2x2_ck *v, matrix_2x2_ck *u_minus_v) {
    //  Calculates the difference between matrix u and matrix v
    //  and returns result in u_minus_v

    for (SU2_color_index a = 0; a < 4; a++) {

        u_minus_v -> m[a] = (u -> m[a]) 
                          - (v -> m[a]);

    }
}

double SU2_trace(const matrix_2x2_ck *u) {
    //	Calculates trace of u

    return 2.0 * (u -> m[0]);
    //	In the Cayley-Klein representation, the trace is
    //	twice the 0th component.
}

double determinant_2x2(const matrix_2x2_ck *u) {
    //  Calculates the determinant of the matrix u

    double det_u = 0.0;

    for (SU2_color_index i = 0; i <= 3; i++) {

        det_u += pow2(u -> m[i]);
        //	In the Cayley-Klein representation, the determinant
        //	is the sum of the squares of the components.

    }

    return det_u;
}

void hermitean_conjugate_2x2(const matrix_2x2_ck *u, matrix_2x2_ck *u_dagger) {
    // Calculates the hermitean conjugate to u
    // and returns result in u_dagger.

    u_dagger -> m[0] = u -> m[0];
    //	In the Cayley-Klein representation, the 0th
    //	component of the conjugate is the same...

    for (SU2_color_index i = 1; i <= 3; i++) {

        u_dagger -> m[i] = -(u -> m[i]);
        //	And the 1, 2 and 3 components are the
        //	same up to a minus sign.

    }
}

void multiplication_by_scalar_2x2(const matrix_2x2_ck *u, const double alpha, matrix_2x2_ck *alpha_times_u) {
    //  Calculates multiplication of Cayley-Klein 2x2 matrix u by scalar alpha
    //  and returns result in alpha_times_u.

    for (SU2_color_index a = 0; a < 4; a++) {

        alpha_times_u -> m[a] = alpha * (u -> m[a]);
        //	Mutiplying each entry.

    }
}

static double SU2_inner_prod(const matrix_2x2_ck *u, const matrix_2x2_ck *v) {
    //	Calculates the "scalar product" of two Cayley-Klein 2x2 matrices
    //	in the Cayley-Klein representation.
    //	Used in the product of two Cayley-Klein 2x2 matrices.

    double inner_prod = (u -> m[0]) * (v -> m[0]);
    //	The 0th component has a plus sign ...

    for (SU2_color_index b = 1; b < 4; b++) {

        inner_prod += -(u -> m[b]) 
                     * (v -> m[b]);

    }
    // and the 1, 2 and 3 components have a minus sign.

    return inner_prod;
}

static void SU2_outer_product(const matrix_2x2_ck *u, const matrix_2x2_ck *v, matrix_2x2_ck *outer_product) {
    //	Calculates the "outer product" of two Cayley-Klein 2x2 matrices
    //	in the Cayley-Klein representation, and returns result in outer_product
    //	Used in the product of two Cayley-Klein 2x2 matrices.

    //	Actually, product taken from the 1, 2 and 3 components only,
    //	in the usual way, as in regular linear algebra.

    outer_product -> m[0] = 0.0;
    outer_product -> m[1] = (u -> m[2]) * (v -> m[3]) 
                          - (u -> m[3]) * (v -> m[2]);
    outer_product -> m[2] = (u -> m[3]) * (v -> m[1]) 
                          - (u -> m[1]) * (v -> m[3]);
    outer_product -> m[3] = (u -> m[1]) * (v -> m[2]) 
                          - (u -> m[2]) * (v -> m[1]);
}

void product_2x2(const matrix_2x2_ck *u, const matrix_2x2_ck *v, matrix_2x2_ck *uv) {
    // Calculates product of 2 Cayley-Klein 2x2 matrices u e v
    // and returns result in uv.

    matrix_2x2_ck u_cross_v;

    uv -> m[0] = SU2_inner_prod(u, v);
    //	In the Cayley-Klein representation, the 0th
    //	components is given by the inner product...

    SU2_outer_product(u, v, &u_cross_v);

    for (SU2_color_index a = 1; a <= 3; a++) {

        uv -> m[a] = (u -> m[a]) * (v -> m[0]) 
                   + (u -> m[0]) * (v -> m[a])
                   - (u_cross_v.m[a]);
        
    }
    //	... and the 1, 2 e 3 components are given
    //	by minus the cross product summed with
    //	a term which mixes 0 and 1, 2, 3 components.
}

void product_three_2x2(const matrix_2x2_ck *u, const matrix_2x2_ck *v, const matrix_2x2_ck *w, matrix_2x2_ck *uvw) {
    //  Calculates product of 3 Cayley-Klein 2x2 matrices u, v and w
    //  and returns result in uvw.

    matrix_2x2_ck uv;
  
    product_2x2(u, v, &uv);
    product_2x2(&uv, w, uvw);

}

void product_four_2x2(const matrix_2x2_ck *u, const matrix_2x2_ck *v, const matrix_2x2_ck *w, const matrix_2x2_ck *x, matrix_2x2_ck *uvwx) {
    //  Calculates product of 4 Cayley-Klein 2x2 matrices u, v, w and x
    //  and returns result in uvwx.

    matrix_2x2_ck uvw;

    product_three_2x2(u, v, w, &uvw);
    product_2x2(&uvw, x, uvwx);

}

void SU2_projection(matrix_2x2_ck *u) {
    //	Projects matrix a to the group SU(2) returning SU(2) matrix a_SU2.

    matrix_2x2_ck u_SU2;

    multiplication_by_scalar_2x2(u, 1.0 / sqrt(determinant_2x2(u)), &u_SU2);
    copy_2x2(&u_SU2, u);

}