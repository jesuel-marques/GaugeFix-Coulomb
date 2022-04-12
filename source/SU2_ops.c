#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>  //	Standard C header files

#include "math_ops.h"  //	Math operations
#include "lattice.h"

// All matrices are in the Cayley-Klein representation
//	u=u[0] SU2_identity + i sum_i=1^3 u[i]sigma[i]
//	where sigma[i] are the Pauli matrices.

void SU2_copy(const double *u, double *u_copy) {
    // Copies u to u_copy

    for (unsigned short a = 0; a < 4; a++) {

        u_copy[a] = u[a];

    }
}

void SU2_set_to_null(double complex *u) {
    // Sets u to be the null matrix in SU(2)

    for (unsigned short a = 0; a < 4; a++) {

        u[a] = 0.0;

    }
}

void SU2_set_to_identity(double complex *u) {
    // Sets u to be the identity matrix in SU(2)

    u[0] = 1.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 0.0;
}

void SU2_accumulate(const double *u, double *acc) {
    // Accumulates the value of u into acc

    for (unsigned short a = 0; a < 4; a++) {

        acc[a] += u[a];

    }
}

void SU2_subtraction(const double *u, const double *v, double *u_minus_v) {
    //  Calculates the difference between matrix u and matrix v
    //  and returns result in u_minus_v

    for (unsigned short a = 0; a < 4; a++) {

        u_minus_v[a] = u[a] - v[a];

    }
}

double SU2_trace(const double *u) {
    //	Calculates trace of u

    return 2.0 * u[0];
    //	In the Cayley-Klein representation, the trace is
    //	twice the 0th component.
}

double SU2_determinant(const double *u) {
    //  Calculates the determinant of the matrix u

    double det_u = 0.0;

    for (unsigned short i = 0; i <= 3; i++) {

        det_u += pow2(u[i]);
        //	In the Cayley-Klein representation, the determinant
        //	is the sum of the squares of the components.

    }

    return det_u;
}

void SU2_hermitean_conjugate(const double *u, double *u_dagger) {
    // Calculates the hermitean conjugate to u
    // and returns result in u_dagger.

    u_dagger[0] = u[0];
    //	In the Cayley-Klein representation, the 0th
    //	component of the conjugate is the same...

    for (unsigned short i = 1; i <= 3; i++) {

        u_dagger[i] = -u[i];
        //	And the 1, 2 and 3 components are the
        //	same up to a minus sign.

    }
}

void SU2_multiplication_by_scalar(const double *u, const double alpha, double *alpha_times_u) {
    //  Calculates multiplicatoin of SU(2) matrix u by scalar alpha
    //  and returns result in alpha_times_u.

    for (unsigned short a = 0; a < 4; a++) {

        alpha_times_u[a] = alpha * u[a];
        //	Mutiplying each entry.

    }
}

static double SU2_inner_prod(const double *u, const double *v) {
    //	Calculates the "scalar product" of two SU(2) matrices
    //	in the Cayley-Klein representation.
    //	Used in the product of two SU(2) matrices.

    double inner_prod = u[0] * v[0];
    //	The 0th component has a plus sign ...

    for (unsigned short b = 1; b < 4; b++) {

        inner_prod += -u[b] * v[b];

    }
    // and the 1, 2 and 3 components have a minus sign.

    return inner_prod;
}

static void SU2_outer_product(const double *u, const double *v, double *outer_product) {
    //	Calculates the "outer product" of two SU(2) matrices
    //	in the Cayley-Klein representation, and returns result in outer_product
    //	Used in the product of two SU(2) matrices.

    //	Actually, product taken from the 1, 2 and 3 components only,
    //	in the usual way, as in regular linear algebra.

    outer_product[0] = 0.0;
    outer_product[1] = u[2] * v[3] - u[3] * v[2];
    outer_product[2] = u[3] * v[1] - u[1] * v[3];
    outer_product[3] = u[1] * v[2] - u[2] * v[1];
}

void SU2_product(const double *u, const double *v, double *uv) {
    // Calculates product of 2 SU(2) matrices u e v
    // and returns result in uv.

    double u_cross_v[4];

    uv[0] = SU2_inner_prod(u, v);
    //	In the Cayley-Klein representation, the 0th
    //	components is given by the inner product...

    SU2_outer_product(u, v, u_cross_v);

    for (unsigned short a = 1; a <= 3; a++) {

        uv[a] = u[a] * v[0] + u[0] * v[a] - *(u_cross_v + a);
        
    }
    //	... and the 1, 2 e 3 components are given
    //	by minus the cross product summed with
    //	a term which mixes 0 and 1, 2, 3 components.
}

void SU2_product_three(const double *u, const double *v, const double *w, double *uvw) {
    //  Calculates product of 3 SU(3) matrices u, v and w
    //  and returns result in uvw.

    double uv[4];
  
    SU2_product(u, v, uv);
    SU2_product(uv, w, uvw);

}

void SU2_product_four(const double *u, const double *v, const double *w, const double *x, double *uvwx) {
    //  Calculates product of 4 SU(3) matrices u, v, w and x
    //  and returns result in uvwx.

    double uvw[4];

    SU2_product_three(u, v, w, uvw);
    SU2_product(uvw, x, uvwx);

}

void SU2_projection(double *u) {
    //	Projects matrix a to the group SU(2) returning SU(2) matrix a_SU2.

    double u_SU2[4];

    SU2_multiplication_by_scalar(u, 1.0 / sqrt(SU2_determinant(u)), u_SU2);
    SU2_copy(u_SU2, u);

}