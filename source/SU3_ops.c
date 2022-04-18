#include <complex.h>
#include <math.h>
#include <stdio.h>  //	Standard header files in C
#include <stdlib.h>

#include "math_ops.h"  //	Math operations
#include "lattice.h"
#include "SU3_ops.h"

void SU3_print_matrix(const double complex *u, const char *name) {
    // Prints the matrix on screen

    printf("\n\n %s \n", name);

    printf("{");
    for (SU3_color_index a = 0; a < Nc; a++) {
        printf("{");

        for (SU3_color_index b = 0; b < Nc; b++) {
            printf("%.16lf+I(%.16lf)", creal(u[Nc * a + b]), cimag(u[Nc * a + b]));
            
            b != 2 ?  printf(",") : 0 ;
        }

        a != 2 ?  printf("},\n") : 0 ;
    }

    printf("}}\n\n");

    getchar();
}

void SU3_copy(const double complex *u, double complex *u_copy) {
    // Copies u to u_copy

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u_copy[a * Nc + b] = u[a * Nc + b];

        }
    }
}

void SU3_convert_fd(const float complex *u_float, double complex *u_double) {
    // Converts link with single precision u_float to u_double with double precision

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u_double[a * Nc + b] = (double complex) u_float[a * Nc + b];

        }
    }
}

void SU3_convert_df(const double complex *u_double, float complex *u_float) {
    // Converts link with single precision u_float to u_double with double precision

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u_float[a * Nc + b] = (float complex) u_double[a * Nc + b];

        }
    }
}

void SU3_set_to_null(double complex *u) {
    // Sets u to be the null matrix in SU(3)

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u[a * Nc + b] = 0.0;

        }
    }
}

void SU3_set_to_identity(double complex *u) {
    // Sets u to be the identity matrix in SU(3)

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

                u[a * Nc + b] = a != b ? 0.0 : 1.0 ;

        }
    }
}

void SU3_accumulate(const double complex *u, double complex *acc) {
    // Accumulates the value of u into acc

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            acc[a * Nc + b] += u[a * Nc + b];
        
        }
    }
}

void SU3_subtraction(const double complex *u, const double complex *v,
                     double complex *u_minus_v) {
    //  Calculates the difference between matrix u and matrix v
    //  and returns result in u_minus_v

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u_minus_v[a * Nc + b] = u[a * Nc + b] - v[a * Nc + b];

        }
    }
}

double complex SU3_trace(const double complex *u) {
    //	Calculates the trace of the matrix u
    //  and returns result (complex) in tr

    return u[0 * Nc + 0] + u[1 * Nc + 1] + u[2 * Nc + 2];
}

double complex SU3_determinant(const double complex *u) {
    //  Calculates the determinant of the matrix u
    //  and returns result, a complex number, in det.

    //  The determinant is calculated in the usual manner
    //  using Leibniz formula

    return u[0 * Nc + 0] * (u[1 * Nc + 1] * u[2 * Nc + 2] - u[1 * Nc + 2] * u[2 * Nc + 1]) 
         + u[0 * Nc + 1] * (u[1 * Nc + 2] * u[2 * Nc + 0] - u[1 * Nc + 0] * u[2 * Nc + 2]) 
         + u[0 * Nc + 2] * (u[1 * Nc + 0] * u[2 * Nc + 1] - u[1 * Nc + 1] * u[2 * Nc + 0]);
}

void SU3_hermitean_conjugate(const double complex *u, double complex *u_dagger) {
    // Calculates the hermitean conjugate to u
    // and returns result in u_dagger.

    u_dagger[0 * Nc + 0] = conj(u[0 * Nc + 0]);
    u_dagger[1 * Nc + 1] = conj(u[1 * Nc + 1]);
    u_dagger[2 * Nc + 2] = conj(u[2 * Nc + 2]);

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < a; b++) {
            if (a != b) {

                u_dagger[b * Nc + a] = conj(u[a * Nc + b]);
                u_dagger[a * Nc + b] = conj(u[b * Nc + a]);

            }
        }
    }
}

void SU3_multiplication_by_scalar(const double complex alpha, const double complex *u,
                                  double complex *alpha_times_u) {
    //  Calculates multiplication of SU(3) matrix u by scalar alpha
    //  and returns result in alphatimesu.

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            alpha_times_u[a * Nc + b] = alpha * u[a * Nc + b];
            //	Mutiplying each entry.

        }
    }
}

void SU3_substitution_multiplication_by_scalar(const double complex alpha, double complex *u) {
    //  Calculates multiplication of SU(3) matrix u by scalar alpha
    //  and returns result in u.

    for (SU3_color_index a = 0; a < Nc; a++) {
        for (SU3_color_index b = 0; b < Nc; b++) {

            u[a * Nc + b] *= alpha;
            //	Mutiplying each entry.
            
        }
    }
}

void SU3_product(const double complex *u, const double complex *v, double complex *uv) {
    // Calculates product of 2 SU(3) matrices u e v
    // and returns result in uv.

    for (SU3_color_index a = 0; a < Nc; a++) {          //  lines
        for (SU3_color_index b = 0; b < Nc; b++) {      //  columns
            uv[a * Nc + b] = 0.0;
            for (SU3_color_index c = 0; c < Nc; c++) {  //  dummy index

                uv[a * Nc + b] += u[a * Nc + c] * v[c * Nc + b];
                //  Usual matrix multiplication.
            }
        }
    }
}

void SU3_product_three(const double complex *u, const double complex *v, const double complex *w,
                       double complex *uvw) {
    //  Calculates product of 3 SU(3) matrices u, v and w
    //  and returns result in uvw.

    for (SU3_color_index a = 0; a < Nc; a++) {              //  lines
        for (SU3_color_index b = 0; b < Nc; b++) {          //  columns
            uvw[a * Nc + b] = 0.0;
            for (SU3_color_index c = 0; c < Nc; c++) {      //  dummy index 1
                for (SU3_color_index e = 0; e < Nc; e++) {  //  dummy index 2

                    uvw[a * Nc + b] += u[a * Nc + c] * v[c * Nc + e] * w[e * Nc + b];
                    //  Usual matrix multiplication.
                }
            }
        }
    }
}

void SU3_product_four(const double complex *u, const double complex *v, const double complex *w, const double complex *x, double complex *uvwx){
    //  Calculates product of 3 SU(3) matrices u, v and w
    //  and returns result in uvw.

    for (SU3_color_index a = 0; a < Nc; a++) {              //  lines
        for (SU3_color_index b = 0; b < Nc; b++) {          //  columns
            uvwx[a * Nc + b] = 0.0;
            for (SU3_color_index c = 0; c < Nc; c++) {            //  dummy index 1
                for (SU3_color_index e = 0; e < Nc; e++) {        //  dummy index 2
                    for (SU3_color_index f = 0; f < Nc; f++) {    //  dummy index 3

                    uvwx[a * Nc + b] += u[a * Nc + c] * v[c * Nc + e] * w[e * Nc + f] * x[f * Nc + b];
                    //  Usual matrix multiplication.
                    }
                }
            }
        }
    }
}

void SU3_accumulate_left_product(const double complex *g, double complex *acc_prod) {
    //	Calculates matrix product between g and acc_prod and accumulates result in acc_prod
    double complex aux_prod[Nc * Nc];
 
    SU3_copy(acc_prod, aux_prod);
    SU3_product(g, aux_prod, acc_prod);

}

void SU3_accumulate_right_product(double complex *acc_prod, const double complex *g) {
    //	Calculates matrix product between acc_prod and g and accumulates result in acc_prod

    double complex aux_prod[Nc * Nc];

    SU3_copy(acc_prod, aux_prod);
    SU3_product(aux_prod, g, acc_prod);

}

void SU3_projection(double complex *x) {
    //	Projects matrix u to the group SU(3) returning SU(3) matrix x at the end.
    //	Follows method found in Gattringer around Eq. 4.27.
    //  More explanation below.

    double sum_absvalue = 0.0;  //  To calculate first two lines of
                                //  projected matrix.

    for (SU3_color_index b = 0; b < Nc; b++) {
        sum_absvalue += pow2(cabsl(x[0 * Nc + b]));
        //	Absvalue of first line of x matrix
        //  to be used for normalization below.
    }

    double complex x_SU3[Nc * Nc];

    //  Used in the Gram-Schmidt
    double complex v_unewconj = 0.0;  //  method to calculate
                                      //  second line of projected
                                      //  matrix.
    
    double complex u_new_conj[Nc];  //  To calculate last line of
    double complex v_new_conj[Nc];  //  projected matrix.

    for (SU3_color_index b = 0; b < Nc; b++) {

        x_SU3[0 * Nc + b] = x[0 * Nc + b] / sqrt(sum_absvalue);
        //	Calculates u_new, first line of x_SU3, projected from x
        //  simply by dividing its first line by its absolute value.

        u_new_conj[b] = conj(x_SU3[0 * Nc + b]);
        //	u_new_conj is u_new's conjugate.

        v_unewconj += x[1 * Nc + b] * u_new_conj[b];
        //	Projection of v, second line of x, onto u_new,
        //  to be used below for the Gram-Schmidt method.

    }

    double complex v_prime[Nc];
    //	Used in the Gram-Schmidt  method to calculate
    //  second line of projected matrix.

    //	Gram-Schmidt method

    for (SU3_color_index b = 0; b < Nc; b++) {

        v_prime[b] = x[1 * Nc + b] - x_SU3[0 * Nc + b] * v_unewconj;
        // Subtraction of the part parallel to u_new of v.

    }

    sum_absvalue = 0.0;

    for (SU3_color_index b = 0; b < Nc; b++) {

        sum_absvalue += pow2(cabsl(v_prime[b]));
        //	Absvalue of second line of projected matrix
        //  before being normalized.

    }

    for (SU3_color_index b = 0; b < Nc; b++) {

        x_SU3[1 * Nc + b] = v_prime[b] / sqrt(sum_absvalue);
        //	Calculates  v_new, second line of the projected matrix
        //  from v_prime, by normalizing it.

        v_new_conj[b] = conj(x_SU3[1 * Nc + b]);
        //	v_new_conj is v_new's conjugate.
        
    }

    //	The third line is calculated by taking the cross product between
    //  u_new_conj and v_new_conj

    //  First component:
    x_SU3[2 * Nc + 0] = (u_new_conj[1] * v_new_conj[2] - u_new_conj[2] * v_new_conj[1]);

    //  Second component:
    x_SU3[2 * Nc + 1] = (u_new_conj[2] * v_new_conj[0] - u_new_conj[0] * v_new_conj[2]);

    //  Third component:
    x_SU3[2 * Nc + 2] = (u_new_conj[0] * v_new_conj[1] - u_new_conj[1] * v_new_conj[0]);

    SU3_copy(x_SU3, x);

}

void SU3_decompose_algebra(const double complex *a, double *a_components) {
    //  Decomposes A in the components of the alfebra of SU(3)
    //  using the Gell-Mann matrices a basis and following the conventions
    //  of Gattringer. The components are then returned in a_components

    // Formulas for the coefficients obtained in Mathematica

    a_components[0] = 0.0;

    a_components[1] = ( creal(a[0 * Nc + 1]) + creal(a[1 * Nc + 0]));

    a_components[2] = (-cimag(a[0 * Nc + 1]) + cimag(a[1 * Nc + 0]));

    a_components[3] = ( creal(a[0 * Nc + 0]) - creal(a[1 * Nc + 1]));

    a_components[4] = ( creal(a[0 * Nc + 2]) + creal(a[2 * Nc + 0]));

    a_components[5] = (-cimag(a[0 * Nc + 2]) + cimag(a[2 * Nc + 0]));

    a_components[6] = ( creal(a[1 * Nc + 2]) + creal(a[2 * Nc + 1]));

    a_components[7] = (-cimag(a[1 * Nc + 2]) + cimag(a[2 * Nc + 1]));

    a_components[8] = ( creal(a[0 * Nc + 0]) + creal(a[1 * Nc + 1]) - 2.0 * creal(a[2 * Nc + 2])) / pow(3.0, 0.5);
}