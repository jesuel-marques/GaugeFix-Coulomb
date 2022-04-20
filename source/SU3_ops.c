#include <complex.h>
#include <math.h>
#include <stdio.h>  //	Standard header files in C
#include <stdlib.h>

#include "math_ops.h"  //	Math operations
#include "lattice.h"
#include "SU3_ops.h"

void print_matrix_3x3(const mtrx_3x3_double *u, const char *name) {
    // Prints the matrix on screen

    printf("\n\n %s \n", name);

    printf("{");
    for (SU3_color_idx a = 0; a < Nc; a++) {
        printf("{");

        for (SU3_color_idx b = 0; b < Nc; b++) {
            printf("%.16lf+I(%.16lf)", creal(u->m[elm(a, b)]), 
                                       cimag(u->m[elm(a, b)]));
            
            b != 2 ?  printf(",") : 0 ;
        }

        a != 2 ?  printf("},\n") : 0 ;
    }

    printf("}}\n\n");

    getchar();
}

void copy_3x3(const mtrx_3x3_double *u, mtrx_3x3_double *u_copy) {
    // Copies u to u_copy

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < Nc; b++) {

            u_copy -> m[elm(a, b)] = u -> m[elm(a, b)];

        }
    }
}

void convert_fd_3x3(const mtrx_3x3_float *u_float, mtrx_3x3_double *u_double) {
    // Converts 3x3 matrix with single precision u_float to u_double with double precision

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < Nc; b++) {

            u_double -> m[elm(a, b)] = (double complex) u_float -> m[elm(a, b)];

        }
    }
}

void convert_df_3x3(const mtrx_3x3_double *u_double, 
                          mtrx_3x3_float *u_float) {
    // Converts 3x3 matrix with single precision 
    // u_float to u_double with double precision

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < Nc; b++) {

            u_float -> m[elm(a, b)] = (float complex) u_double -> m[elm(a, b)];

        }
    }
}

inline void set_null_3x3(mtrx_3x3_double * u) {
    // Sets u to be the 3x3 null matrix

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < Nc; b++) {

            u -> m[elm(a, b)] = 0.0;

        }
    }
}

inline void set_identity_3x3(mtrx_3x3_double *u) {
    // Sets u to be the identity matrix in SU(3)

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < Nc; b++) {

                u -> m[elm(a, b)] = a != b ? 0.0 : 1.0 ;

        }
    }
}

void accumulate_3x3(const mtrx_3x3_double *u, mtrx_3x3_double *acc) {
    // Accumulates the value of u into acc

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < Nc; b++) {

            acc -> m[elm(a, b)] += u -> m[elm(a, b)];
        
        }
    }
}

void subtraction_3x3(const mtrx_3x3_double *u, 
                     const mtrx_3x3_double *v,
                                mtrx_3x3_double *u_minus_v) {
    //  Calculates the difference between matrix u and matrix v
    //  and returns result in u_minus_v

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < Nc; b++) {

            u_minus_v -> m[elm(a, b)] = u -> m[elm(a, b)] 
                                      - v -> m[elm(a, b)];

        }
    }
}

double complex trace_3x3(const mtrx_3x3_double *u) {
    //	Calculates the trace of the matrix u
    //  and returns result (complex) in tr

    return u -> m[elm(0, 0)] 
         + u -> m[elm(1, 1)] 
         + u -> m[elm(2, 2)];
}

inline double complex determinant_3x3(const mtrx_3x3_double *u) {
    //  Calculates the determinant of the matrix u
    //  and returns result, a complex number, in det.

    //  The determinant is calculated in the usual manner
    //  using Leibniz formula

    return u -> m[elm(0, 0)] * (u -> m[elm(1, 1)] * u -> m[elm(2, 2)] 
                              - u -> m[elm(1, 2)] * u -> m[elm(2, 1)]) 
         + u -> m[elm(0, 1)] * (u -> m[elm(1, 2)] * u -> m[elm(2, 0)] 
                              - u -> m[elm(1, 0)] * u -> m[elm(2, 2)]) 
         + u -> m[elm(0, 2)] * (u -> m[elm(1, 0)] * u -> m[elm(2, 1)] 
                              - u -> m[elm(1, 1)] * u -> m[elm(2, 0)]);
}

inline void SU3_hermitean_conjugate(const mtrx_3x3_double *u, mtrx_3x3_double *u_dagger) {
    // Calculates the hermitean conjugate to u
    // and returns result in u_dagger.

    u_dagger -> m[elm(0, 0)] = conj(u -> m[elm(0, 0)]);
    u_dagger -> m[elm(1, 1)] = conj(u -> m[elm(1, 1)]);
    u_dagger -> m[elm(2, 2)] = conj(u -> m[elm(2, 2)]);

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < a; b++) {
            if (a != b) {

                u_dagger -> m[elm(b, a)] = conj(u -> m[elm(a, b)]);
                u_dagger -> m[elm(a, b)] = conj(u -> m[elm(b, a)]);

            }
        }
    }
}

inline void multiplication_by_scalar_3x3(const double complex alpha, const mtrx_3x3_double *u,
                                            mtrx_3x3_double *alpha_times_u) {
    //  Calculates multiplication of 3x3 matrix u by scalar alpha
    //  and returns result in alphatimesu.

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < Nc; b++) {

            alpha_times_u -> m[elm(a, b)] = alpha * u -> m[elm(a, b)];
            //	Mutiplying each entry.

        }
    }
}

inline void substitution_multiplication_by_scalar_3x3(const double complex alpha, 
                                                                mtrx_3x3_double *u) {
    //  Calculates multiplication of 3x3 matrix u by scalar alpha
    //  and returns result in u.

    for (SU3_color_idx a = 0; a < Nc; a++) {
        for (SU3_color_idx b = 0; b < Nc; b++) {

            u -> m[elm(a, b)] *= alpha;
            //	Mutiplying each entry.
            
        }
    }
}

inline void product_3x3(const mtrx_3x3_double *u, 
                        const mtrx_3x3_double *v, 
                                    mtrx_3x3_double *uv) {
    // Calculates product of 2 3x3 matrices u e v
    // and returns result in uv.

    for (SU3_color_idx a = 0; a < Nc; a++) {          //  lines
        for (SU3_color_idx b = 0; b < Nc; b++) {      //  columns
            uv -> m[elm(a, b)] = 0.0;
            for (SU3_color_idx c = 0; c < Nc; c++) {  //  dummy index

                uv -> m[elm(a, b)] += u -> m[elm(a, c)] 
                                    * v -> m[elm(c, b)];
                //  Usual matrix multiplication.
            }
        }
    }
}

void product_three_3x3(const mtrx_3x3_double *u, 
                       const mtrx_3x3_double *v,
                       const mtrx_3x3_double *w,
                                mtrx_3x3_double *uvw) {
    //  Calculates product of 3 3x3 matrices u, v and w
    //  and returns result in uvw.

    for (SU3_color_idx a = 0; a < Nc; a++) {              //  lines
        for (SU3_color_idx b = 0; b < Nc; b++) {          //  columns
            uvw -> m[elm(a, b)] = 0.0;
            for (SU3_color_idx c = 0; c < Nc; c++) {      //  dummy index 1
                for (SU3_color_idx e = 0; e < Nc; e++) {  //  dummy index 2

                    uvw -> m[elm(a, b)] += (u -> m[elm(a, c)]) 
                                         * (v -> m[elm(c, e)]) 
                                         * (w -> m[elm(e, b)]);
                    //  Usual matrix multiplication.
                }
            }
        }
    }
}

void product_four_3x3(const mtrx_3x3_double *u, 
                      const mtrx_3x3_double *v,
                      const mtrx_3x3_double *w,
                      const mtrx_3x3_double *x, 
                                mtrx_3x3_double *uvwx){
    //  Calculates product of 3 3x3 matrices u, v and w
    //  and returns result in uvw.

    for (SU3_color_idx a = 0; a < Nc; a++) {              //  lines
        for (SU3_color_idx b = 0; b < Nc; b++) {          //  columns
            uvwx -> m[elm(a, b)] = 0.0;
            for (SU3_color_idx c = 0; c < Nc; c++) {            //  dummy index 1
                for (SU3_color_idx e = 0; e < Nc; e++) {        //  dummy index 2
                    for (SU3_color_idx f = 0; f < Nc; f++) {    //  dummy index 3

                    uvwx -> m[elm(a, b)] += (u -> m[elm(a, c)]) 
                                          * (v -> m[elm(c, e)]) 
                                          * (w -> m[elm(e, f)]) 
                                          * (x -> m[elm(f, b)]);
                    //  Usual matrix multiplication.
                    }
                }
            }
        }
    }
}

inline void accumulate_left_prod_3x3(const mtrx_3x3_double *g, mtrx_3x3_double *acc_prod) {
    //	Calculates matrix product between g and acc_prod and accumulates result in acc_prod
    mtrx_3x3_double aux_prod;
 
    copy_3x3(acc_prod, &aux_prod);
    product_3x3(g, &aux_prod, acc_prod);

}

inline void accumulate_right_prod_3x3(mtrx_3x3_double *acc_prod, const mtrx_3x3_double *g) {
    //	Calculates matrix product between acc_prod and g and accumulates result in acc_prod

    mtrx_3x3_double aux_prod;

    copy_3x3(acc_prod, &aux_prod);
    product_3x3(&aux_prod, g, acc_prod);

}

inline void projection_SU3(mtrx_3x3_double *x) {
    //	Projects matrix u to the group SU(3) returning SU(3) matrix x at the end.
    //	Follows method found in Gattringer around Eq. 4.27.
    //  More explanation below.

    double sum_absvalue = 0.0;  //  To calculate first two lines of
                                //  projected matrix.

    for (SU3_color_idx b = 0; b < Nc; b++) {
        sum_absvalue += pow2(cabsl(x -> m[elm(0, b)]));
        //	Absvalue of first line of x matrix
        //  to be used for normalization below.
    }

    mtrx_3x3_double x_SU3;

    //  Used in the Gram-Schmidt
    double complex v_unewconj = 0.0;  //  method to calculate
                                      //  second line of projected
                                      //  matrix.
    
    color_3_vec u_new_conj;  //  To calculate last line of
    color_3_vec v_new_conj;  //  projected matrix.

    for (SU3_color_idx b = 0; b < Nc; b++) {

        x_SU3.m[elm(0, b)] = (x -> m[elm(0, b)]) / sqrt(sum_absvalue);
        //	Calculates u_new, first line of x_SU3, projected from x
        //  simply by dividing its first line by its absolute value.

        u_new_conj.m[b] = conj(x_SU3.m[elm(0, b)]);
        //	u_new_conj is u_new's conjugate.

        v_unewconj += (x -> m[elm(1, b)]) * u_new_conj.m[b];
        //	Projection of v, second line of x, onto u_new,
        //  to be used below for the Gram-Schmidt method.

    }

    double complex v_prime[Nc];
    //	Used in the Gram-Schmidt  method to calculate
    //  second line of projected matrix.

    //	Gram-Schmidt method

    for (SU3_color_idx b = 0; b < Nc; b++) {

        v_prime[b] = (x -> m[elm(1, b)]) - x_SU3.m[elm(0, b)] * v_unewconj;
        // Subtraction of the part parallel to u_new of v.

    }

    sum_absvalue = 0.0;

    for (SU3_color_idx b = 0; b < Nc; b++) {

        sum_absvalue += pow2(cabsl(v_prime[b]));
        //	Absvalue of second line of projected matrix
        //  before being normalized.

    }

    for (SU3_color_idx b = 0; b < Nc; b++) {

        x_SU3.m[elm(1, b)] = v_prime[b] / sqrt(sum_absvalue);
        //	Calculates  v_new, second line of the projected matrix
        //  from v_prime, by normalizing it.

        v_new_conj.m[b] = conj(x_SU3.m[elm(1, b)]);
        //	v_new_conj is v_new's conjugate.
        
    }

    //	The third line is calculated by taking the cross product between
    //  u_new_conj and v_new_conj

    //  First component:
    x_SU3.m[elm(2, 0)] = (u_new_conj.m[1] * v_new_conj.m[2] 
                        - u_new_conj.m[2] * v_new_conj.m[1]);

    //  Second component:
    x_SU3.m[elm(2, 1)] = (u_new_conj.m[2] * v_new_conj.m[0] 
                        - u_new_conj.m[0] * v_new_conj.m[2]);

    //  Third component:
    x_SU3.m[elm(2, 2)] = (u_new_conj.m[0] * v_new_conj.m[1] 
                        - u_new_conj.m[1] * v_new_conj.m[0]);

    copy_3x3(&x_SU3, x);

}

void decompose_algebra_SU3(const mtrx_3x3_double *a, matrix_SU3_alg *a_components) {
    //  Decomposes A in the components of the alfebra of SU(3)
    //  using the Gell-Mann matrices a basis and following the conventions
    //  of Gattringer. The components are then returned in a_components

    // Formulas for the coefficients obtained in Mathematica

    a_components -> m[0] = 0.0;

    a_components -> m[1] = ( creal(a -> m[elm(0, 1)]) 
                           + creal(a -> m[elm(1, 0)]));

    a_components -> m[2] = (-cimag(a -> m[elm(0, 1)]) 
                           + cimag(a -> m[elm(1, 0)]));

    a_components -> m[3] = ( creal(a -> m[elm(0, 0)]) 
                           - creal(a -> m[elm(1, 1)]));

    a_components -> m[4] = ( creal(a -> m[elm(0, 2)]) 
                           + creal(a -> m[elm(2, 0)]));

    a_components -> m[5] = (-cimag(a -> m[elm(0, 2)]) 
                           + cimag(a -> m[elm(2, 0)]));

    a_components -> m[6] = ( creal(a -> m[elm(1, 2)]) 
                           + creal(a -> m[elm(2, 1)]));

    a_components -> m[7] = (-cimag(a -> m[elm(1, 2)]) 
                           + cimag(a -> m[elm(2, 1)]));

    a_components -> m[8] = ( creal(a -> m[elm(0, 0)]) 
                           + creal(a -> m[elm(1, 1)]) 
                     - 2.0 * creal(a -> m[elm(2, 2)])) / pow(3.0, 0.5);
}