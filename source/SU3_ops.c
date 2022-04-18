#include <complex.h>
#include <math.h>
#include <stdio.h>  //	Standard header files in C
#include <stdlib.h>

#include "math_ops.h"  //	Math operations
#include "lattice.h"
#include "SU3_ops.h"

void SU3_print_matrix(const matrix_3x3_double *u, const char *name) {
    // Prints the matrix on screen

    printf("\n\n %s \n", name);

    printf("{");
    for (SU3_color_index a = 0; a < Nc; a++) {
        printf("{");

        for (SU3_color_index b = 0; b < Nc; b++) {
            printf("%.16lf+I(%.16lf)", creal(u->m[Nc * a + b]), cimag(u->m[Nc * a + b]));
            
            b != 2 ?  printf(",") : 0 ;
        }

        a != 2 ?  printf("},\n") : 0 ;
    }

    printf("}}\n\n");

    getchar();
}

void SU3_copy(const matrix_3x3_double *u, matrix_3x3_double *u_copy) {
    // Copies u to u_copy

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u_copy -> m[a * Nc + b] = u -> m[a * Nc + b];

        }
    }
}

void SU3_convert_fd(const matrix_3x3_float *u_float, matrix_3x3_double *u_double) {
    // Converts link with single precision u_float to u_double with double precision

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u_double -> m[a * Nc + b] = (double complex) u_float -> m[a * Nc + b];

        }
    }
}

void SU3_convert_df(const matrix_3x3_double *u_double, matrix_3x3_float *u_float) {
    // Converts link with single precision u_float to u_double with double precision

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u_float -> m[a * Nc + b] = (float complex) u_double -> m[a * Nc + b];

        }
    }
}

void SU3_set_to_null(matrix_3x3_double * u) {
    // Sets u to be the null matrix in SU(3)

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u -> m[a * Nc + b] = 0.0;

        }
    }
}

void SU3_set_to_identity(matrix_3x3_double *u) {
    // Sets u to be the identity matrix in SU(3)

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

                u -> m[a * Nc + b] = a != b ? 0.0 : 1.0 ;

        }
    }
}

void SU3_accumulate(const matrix_3x3_double *u, matrix_3x3_double *acc) {
    // Accumulates the value of u into acc

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            acc -> m[a * Nc + b] += u -> m[a * Nc + b];
        
        }
    }
}

void SU3_subtraction(const matrix_3x3_double *u, const matrix_3x3_double *v,
                     matrix_3x3_double *u_minus_v) {
    //  Calculates the difference between matrix u and matrix v
    //  and returns result in u_minus_v

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u_minus_v -> m[a * Nc + b] = u -> m[a * Nc + b] - v -> m[a * Nc + b];

        }
    }
}

double complex SU3_trace(const matrix_3x3_double *u) {
    //	Calculates the trace of the matrix u
    //  and returns result (complex) in tr

    return u -> m[0 * Nc + 0] 
         + u -> m[1 * Nc + 1] 
         + u -> m[2 * Nc + 2];
}

double complex SU3_determinant(const matrix_3x3_double *u) {
    //  Calculates the determinant of the matrix u
    //  and returns result, a complex number, in det.

    //  The determinant is calculated in the usual manner
    //  using Leibniz formula

    return u -> m[0 * Nc + 0] * (u -> m[1 * Nc + 1] * u -> m[2 * Nc + 2] - u -> m[1 * Nc + 2] * u -> m[2 * Nc + 1]) 
         + u -> m[0 * Nc + 1] * (u -> m[1 * Nc + 2] * u -> m[2 * Nc + 0] - u -> m[1 * Nc + 0] * u -> m[2 * Nc + 2]) 
         + u -> m[0 * Nc + 2] * (u -> m[1 * Nc + 0] * u -> m[2 * Nc + 1] - u -> m[1 * Nc + 1] * u -> m[2 * Nc + 0]);
}

void SU3_hermitean_conjugate(const matrix_3x3_double *u, matrix_3x3_double *u_dagger) {
    // Calculates the hermitean conjugate to u
    // and returns result in u_dagger.

    u_dagger -> m[0 * Nc + 0] = conj(u -> m[0 * Nc + 0]);
    u_dagger -> m[1 * Nc + 1] = conj(u -> m[1 * Nc + 1]);
    u_dagger -> m[2 * Nc + 2] = conj(u -> m[2 * Nc + 2]);

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < a; b++) {
            if (a != b) {

                u_dagger -> m[b * Nc + a] = conj(u -> m[a * Nc + b]);
                u_dagger -> m[a * Nc + b] = conj(u -> m[b * Nc + a]);

            }
        }
    }

    return u_dagger;
}

void SU3_multiplication_by_scalar(const double complex alpha, const matrix_3x3_double *u,
                                  matrix_3x3_double *alpha_times_u) {
    //  Calculates multiplication of SU(3) matrix u by scalar alpha
    //  and returns result in alphatimesu.

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            alpha_times_u -> m[a * Nc + b] = alpha * u -> m[a * Nc + b];
            //	Mutiplying each entry.

        }
    }
}

void SU3_substitution_multiplication_by_scalar(const double complex alpha, matrix_3x3_double *u) {
    //  Calculates multiplication of SU(3) matrix u by scalar alpha
    //  and returns result in u.

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u -> m[a * Nc + b] *= alpha;
            //	Mutiplying each entry.
            
        }
    }
}

void SU3_product(const matrix_3x3_double *u, const matrix_3x3_double *v, matrix_3x3_double *uv) {
    // Calculates product of 2 SU(3) matrices u e v
    // and returns result in uv.

    for (SU3_color_index a = 0; a < Nc; a++) {          //  lines
        for (SU3_color_index b = 0; b < Nc; b++) {      //  columns
            uv -> m[a * Nc + b] = 0.0;
            for (SU3_color_index c = 0; c < Nc; c++) {  //  dummy index

                uv -> m[a * Nc + b] += u -> m[a * Nc + c] * v -> m[c * Nc + b];
                //  Usual matrix multiplication.
            }
        }
    }
}

void SU3_product_three(const matrix_3x3_double *u, const matrix_3x3_double *v, const matrix_3x3_double *w,
                       matrix_3x3_double *uvw) {
    //  Calculates product of 3 SU(3) matrices u, v and w
    //  and returns result in uvw.

    for (SU3_color_index a = 0; a < Nc; a++) {              //  lines
        for (SU3_color_index b = 0; b < Nc; b++) {          //  columns
            uvw -> m[a * Nc + b] = 0.0;
            for (SU3_color_index c = 0; c < Nc; c++) {      //  dummy index 1
                for (SU3_color_index e = 0; e < Nc; e++) {  //  dummy index 2

                    uvw -> m[a * Nc + b] += u -> m[a * Nc + c] * v -> m[c * Nc + e] * w -> m[e * Nc + b];
                    //  Usual matrix multiplication.
                }
            }
        }
    }
}

void SU3_product_four(const matrix_3x3_double *u, const matrix_3x3_double *v, const matrix_3x3_double *w, const matrix_3x3_double *x, matrix_3x3_double *uvwx){
    //  Calculates product of 3 SU(3) matrices u, v and w
    //  and returns result in uvw.

    for (SU3_color_index a = 0; a < Nc; a++) {              //  lines
        for (SU3_color_index b = 0; b < Nc; b++) {          //  columns
            uvwx -> m[a * Nc + b] = 0.0;
            for (SU3_color_index c = 0; c < Nc; c++) {            //  dummy index 1
                for (SU3_color_index e = 0; e < Nc; e++) {        //  dummy index 2
                    for (SU3_color_index f = 0; f < Nc; f++) {    //  dummy index 3

                    uvwx -> m[a * Nc + b] += u -> m[a * Nc + c] * v -> m[c * Nc + e] * w -> m[e * Nc + f] * x -> m[f * Nc + b];
                    //  Usual matrix multiplication.
                    }
                }
            }
        }
    }
}

void SU3_accumulate_left_product(const matrix_3x3_double *g, matrix_3x3_double *acc_prod) {
    //	Calculates matrix product between g and acc_prod and accumulates result in acc_prod
    matrix_3x3_double aux_prod;
 
    SU3_copy(acc_prod, &aux_prod);
    SU3_product(g, &aux_prod, acc_prod);

}

void SU3_accumulate_right_product(matrix_3x3_double *acc_prod, const matrix_3x3_double *g) {
    //	Calculates matrix product between acc_prod and g and accumulates result in acc_prod

    matrix_3x3_double aux_prod;

    SU3_copy(acc_prod, &aux_prod);
    SU3_product(&aux_prod, g, acc_prod);

}

void SU3_projection(matrix_3x3_double *x) {
    //	Projects matrix u to the group SU(3) returning SU(3) matrix x at the end.
    //	Follows method found in Gattringer around Eq. 4.27.
    //  More explanation below.

    double sum_absvalue = 0.0;  //  To calculate first two lines of
                                //  projected matrix.

    for (SU3_color_index b = 0; b < Nc; b++) {
        sum_absvalue += pow2(cabsl(x -> m[0 * Nc + b]));
        //	Absvalue of first line of x matrix
        //  to be used for normalization below.
    }

    matrix_3x3_double x_SU3;

    //  Used in the Gram-Schmidt
    double complex v_unewconj = 0.0;  //  method to calculate
                                      //  second line of projected
                                      //  matrix.
    
    color_3_vec u_new_conj;  //  To calculate last line of
    color_3_vec v_new_conj;  //  projected matrix.

    for (SU3_color_index b = 0; b < Nc; b++) {

        x_SU3.m[0 * Nc + b] = (x -> m[0 * Nc + b]) / sqrt(sum_absvalue);
        //	Calculates u_new, first line of x_SU3, projected from x
        //  simply by dividing its first line by its absolute value.

        u_new_conj.m[b] = conj(x_SU3.m[0 * Nc + b]);
        //	u_new_conj is u_new's conjugate.

        v_unewconj += (x -> m[1 * Nc + b]) * u_new_conj.m[b];
        //	Projection of v, second line of x, onto u_new,
        //  to be used below for the Gram-Schmidt method.

    }

    double complex v_prime[Nc];
    //	Used in the Gram-Schmidt  method to calculate
    //  second line of projected matrix.

    //	Gram-Schmidt method

    for (SU3_color_index b = 0; b < Nc; b++) {

        v_prime[b] = (x -> m[1 * Nc + b]) - x_SU3.m[0 * Nc + b] * v_unewconj;
        // Subtraction of the part parallel to u_new of v.

    }

    sum_absvalue = 0.0;

    for (SU3_color_index b = 0; b < Nc; b++) {

        sum_absvalue += pow2(cabsl(v_prime[b]));
        //	Absvalue of second line of projected matrix
        //  before being normalized.

    }

    for (SU3_color_index b = 0; b < Nc; b++) {

        x_SU3.m[1 * Nc + b] = v_prime[b] / sqrt(sum_absvalue);
        //	Calculates  v_new, second line of the projected matrix
        //  from v_prime, by normalizing it.

        v_new_conj.m[b] = conj(x_SU3.m[1 * Nc + b]);
        //	v_new_conj is v_new's conjugate.
        
    }

    //	The third line is calculated by taking the cross product between
    //  u_new_conj and v_new_conj

    //  First component:
    x_SU3.m[2 * Nc + 0] = (u_new_conj.m[1] * v_new_conj.m[2] - u_new_conj.m[2] * v_new_conj.m[1]);

    //  Second component:
    x_SU3.m[2 * Nc + 1] = (u_new_conj.m[2] * v_new_conj.m[0] - u_new_conj.m[0] * v_new_conj.m[2]);

    //  Third component:
    x_SU3.m[2 * Nc + 2] = (u_new_conj.m[0] * v_new_conj.m[1] - u_new_conj.m[1] * v_new_conj.m[0]);

    SU3_copy(&x_SU3, x);

}

void SU3_decompose_algebra(const matrix_3x3_double *a, matrix_SU3_alg *a_components) {
    //  Decomposes A in the components of the alfebra of SU(3)
    //  using the Gell-Mann matrices a basis and following the conventions
    //  of Gattringer. The components are then returned in a_components

    // Formulas for the coefficients obtained in Mathematica

    a_components -> m[0] = 0.0;

    a_components -> m[1] = ( creal(a -> m[0 * Nc + 1]) + creal(a -> m[1 * Nc + 0]));

    a_components -> m[2] = (-cimag(a -> m[0 * Nc + 1]) + cimag(a -> m[1 * Nc + 0]));

    a_components -> m[3] = ( creal(a -> m[0 * Nc + 0]) - creal(a -> m[1 * Nc + 1]));

    a_components -> m[4] = ( creal(a -> m[0 * Nc + 2]) + creal(a -> m[2 * Nc + 0]));

    a_components -> m[5] = (-cimag(a -> m[0 * Nc + 2]) + cimag(a -> m[2 * Nc + 0]));

    a_components -> m[6] = ( creal(a -> m[1 * Nc + 2]) + creal(a -> m[2 * Nc + 1]));

    a_components -> m[7] = (-cimag(a -> m[1 * Nc + 2]) + cimag(a -> m[2 * Nc + 1]));

    a_components -> m[8] = ( creal(a -> m[0 * Nc + 0]) + creal(a -> m[1 * Nc + 1]) - 2.0 * creal(a -> m[2 * Nc + 2])) / pow(3.0, 0.5);
}