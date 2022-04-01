#include <complex.h>
#include <math.h>
#include <stdio.h>  //	Standard header files in C
#include <stdlib.h>

#include "math_ops.h"  //	Math operations
#include "lattice.h"

void SU3_print_matrix(const double complex *u, const char *name) {
    // Prints the matrix on screen

    printf("\n\n %s \n", name);

    printf("{");
    for (unsigned short i = 0; i < 3; i++) {
        printf("{");

        for (unsigned short j = 0; j < 3; j++) {
            printf("%.16lf+I(%.16lf)", creal(u[3 * i + j]), cimag(u[3 * i + j]));
            if (j != 2) printf(",");
        }
        if (i != 2) {
            printf("},\n");
        }
    }

    printf("}}");

    getchar();
}

void SU3_copy(const double complex *u, double complex *u_copy) {
    // Copies u to u_copy

    for (unsigned short a = 0; a < 3; a++) {
        for (unsigned short b = 0; b < 3; b++) {
            u_copy[a * 3 + b] = u[a * 3 + b];
        }
    }
}

void SU3_set_to_null(double complex *u) {
    // Sets u to be the null matrix in SU(3)

    for (unsigned short a = 0; a < 3; a++) {
        for (unsigned short b = 0; b < 3; b++) {
            u[a * 3 + b] = 0.0;
        }
    }
}

void SU3_set_to_identity(double complex *u) {
    // Sets u to be the identity matrix in SU(3)

    for (unsigned short a = 0; a < 3; a++) {
        for (unsigned short b = 0; b < 3; b++) {
            if (a != b) {
                u[a * 3 + b] = 0.0;
            } else {
                u[a * 3 + b] = 1.0;
            }
        }
    }
}

void SU3_accumulate(const double complex *u, double complex *acc) {
    // Accumulates the value of u into acc

    for (unsigned short a = 0; a < 3; a++) {
        for (unsigned short b = 0; b < 3; b++) {
            acc[a * 3 + b] += u[a * 3 + b];
        }
    }
}

void SU3_subtraction(const double complex *u, const double complex *v,
                     double complex *u_minus_v) {
    //  Calculates the difference between matrix u and matrix v
    //  and returns result in u_minus_v

    for (unsigned short a = 0; a < 3; a++) {
        for (unsigned short b = 0; b < 3; b++) {
            u_minus_v[a * 3 + b] = u[a * 3 + b] - v[a * 3 + b];
        }
    }
}

double complex SU3_trace(const double complex *u) {
    //	Calculates the trace of the matrix u
    //  and returns result (complex) in tr

    return u[0 * 3 + 0] + u[1 * 3 + 1] + u[2 * 3 + 2];
}

double complex SU3_determinant(const double complex *u) {
    //  Calculates the determinant of the matrix u
    //  and returns result, a complex number, in det.

    //  The determinant is calculated in the usual manner
    //  using Leibniz formula

    return u[0 * 3 + 0] * (u[1 * 3 + 1] * u[2 * 3 + 2] - u[1 * 3 + 2] * u[2 * 3 + 1]) + u[0 * 3 + 1] * (u[1 * 3 + 2] * u[2 * 3 + 0] - u[1 * 3 + 0] * u[2 * 3 + 2]) + u[0 * 3 + 2] * (u[1 * 3 + 0] * u[2 * 3 + 1] - u[1 * 3 + 1] * u[2 * 3 + 0]);
}

void SU3_hermitean_conjugate(const double complex *u, double complex *u_dagger) {
    // Calculates the hermitean conjugate to u
    // and returns result in u_dagger.

    u_dagger[0 * 3 + 0] = conj(u[0 * 3 + 0]);
    u_dagger[1 * 3 + 1] = conj(u[1 * 3 + 1]);
    u_dagger[2 * 3 + 2] = conj(u[2 * 3 + 2]);

    for (unsigned short a = 0; a < 3; a++) {
        for (unsigned short b = 0; b < a; b++) {
            if (a != b) {
                u_dagger[b * 3 + a] = conj(u[a * 3 + b]);
                u_dagger[a * 3 + b] = conj(u[b * 3 + a]);
            }
        }
    }
}

void SU3_multiplication_by_scalar(const double complex alpha, const double complex *u,
                                  double complex *alpha_times_u) {
    //  Calculates multiplication of SU(3) matrix u by scalar alpha
    //  and returns result in alphatimesu.

    for (unsigned short a = 0; a < 3; a++) {
        for (unsigned short b = 0; b < 3; b++) {
            alpha_times_u[a * 3 + b] = alpha * u[a * 3 + b];
            //	Mutiplying each entry.
        }
    }
}

void SU3_substitution_multiplication_by_scalar(const double complex alpha, double complex *u) {
    //  Calculates multiplication of SU(3) matrix u by scalar alpha
    //  and returns result in u.

    for (unsigned short a = 0; a < 3; a++) {
        for (unsigned short b = 0; b < 3; b++) {
            u[a * 3 + b] *= alpha;
            //	Mutiplying each entry.
        }
    }
}

void SU3_product(const double complex *u, const double complex *v, double complex *uv) {
    // Calculates product of 2 SU(3) matrices u e v
    // and returns result in uv.

    SU3_set_to_null(uv);
    for (unsigned short a = 0; a < 3; a++) {          //  lines
        for (unsigned short b = 0; b < 3; b++) {      //  columns
            for (unsigned short c = 0; c < 3; c++) {  //  dummy index

                uv[a * 3 + b] += u[a * 3 + c] * v[c * 3 + b];
                //  Usual matrix multiplication.
            }
        }
    }
}

void SU3_product_three(const double complex *u, const double complex *v, const double complex *w,
                       double complex *uvw) {
    //  Calculates product of 3 SU(3) matrices u, v and w
    //  and returns result in uvw.

    SU3_set_to_null(uvw);
    for (unsigned short a = 0; a < 3; a++) {              //  lines
        for (unsigned short b = 0; b < 3; b++) {          //  columns
            for (unsigned short c = 0; c < 3; c++) {      //  dummy index 1
                for (unsigned short e = 0; e < 3; e++) {  //  dummy index 2

                    uvw[a * 3 + b] += u[a * 3 + c] * v[c * 3 + e] * w[e * 3 + b];
                    //  Usual matrix multiplication.
                }
            }
        }
    }
}

void SU3_accumulate_left_product(const double complex *g, double complex *acc_prod) {
    //	Calculates matrix product between g and acc_prod and accumulates result in acc_prod
    double complex aux_prod[3 * 3];
 
    SU3_copy(acc_prod, aux_prod);
    SU3_product(g, aux_prod, acc_prod);

}

void SU3_accumulate_right_product(double complex *acc_prod, const double complex *g) {
    //	Calculates matrix product between acc_prod and g and accumulates result in acc_prod

    double complex aux_prod[3 * 3];

    SU3_copy(acc_prod, aux_prod);
    SU3_product(aux_prod, g, acc_prod);

}

void SU3_projection(double complex *x) {
    //	Projects matrix u to the group SU(3) returning SU(3) matrix x at the end.
    //	Follows method found in Gattringer around Eq. 4.27.
    //  More explanation below.

    double sum_absvalue = 0.0;  //  To calculate first two lines of
                                //  projected matrix.

    for (unsigned short b = 0; b < 3; b++) {
        sum_absvalue += pow2(cabsl(x[0 * 3 + b]));
        //	Absvalue of first line of x matrix
        //  to be used for normalization below.
    }

    double complex x_SU3[3 * 3];

    //  Used in the Gram-Schmidt
    double complex v_unewconj = 0.0;  //  method to calculate
                                      //  second line of projected
                                      //  matrix.
    
    double complex u_new_conj[3];  //  To calculate last line of
    double complex v_new_conj[3];  //  projected matrix.

    for (unsigned short b = 0; b < 3; b++) {
        x_SU3[0 * 3 + b] = x[0 * 3 + b] / sqrt(sum_absvalue);
        //	Calculates u_new, first line of x_SU3, projected from x
        //  simply by dividing its first line by its absolute value.

        u_new_conj[b] = conj(x_SU3[0 * 3 + b]);
        //	u_new_conj is u_new's conjugate.

        v_unewconj += x[1 * 3 + b] * u_new_conj[b];
        //	Projection of v, second line of x, onto u_new,
        //  to be used below for the Gram-Schmidt method.
    }

    double complex v_prime[3];
    //	Used in the Gram-Schmidt  method to calculate
    //  second line of projected matrix.

    //	Gram-Schmidt method

    for (unsigned short b = 0; b < 3; b++) {
        v_prime[b] = x[1 * 3 + b] - x_SU3[0 * 3 + b] * v_unewconj;
        // Subtraction of the part parallel to u_new of v.
    }

    sum_absvalue = 0.0;

    for (unsigned short b = 0; b < 3; b++) {
        sum_absvalue += pow2(cabsl(v_prime[b]));
        //	Absvalue of second line of projected matrix
        //  before being normalized.
    }

    for (unsigned short b = 0; b < 3; b++) {
        x_SU3[1 * 3 + b] = v_prime[b] / sqrt(sum_absvalue);
        //	Calculates  v_new, second line of the projected matrix
        //  from v_prime, by normalizing it.

        v_new_conj[b] = conj(x_SU3[1 * 3 + b]);
        //	v_new_conj is v_new's conjugate.
    }

    //	The third line is calculated by taking the cross product between
    //  u_new_conj and v_new_conj

    //  First component:
    x_SU3[2 * 3 + 0] = (u_new_conj[1] * v_new_conj[2] - u_new_conj[2] * v_new_conj[1]);

    //  Second component:
    x_SU3[2 * 3 + 1] = (u_new_conj[2] * v_new_conj[0] - u_new_conj[0] * v_new_conj[2]);

    //  Third component:
    x_SU3[2 * 3 + 2] = (u_new_conj[0] * v_new_conj[1] - u_new_conj[1] * v_new_conj[0]);

    SU3_copy(x_SU3, x);

}

void SU3_decompose_algebra(const double complex *a, double *a_components) {
    //  Decomposes A in the components of the alfebra of SU(3)
    //  using the Gell-Mann matrices a basis and following the conventions
    //  of Gattringer. The components are then returned in a_components

    // Formulas for the coefficients obtained in Mathematica

    a_components[0] = 0.0;

    a_components[1] = (creal(a[0 * 3 + 1]) + creal(a[1 * 3 + 0]));

    a_components[2] = (-cimag(a[0 * 3 + 1]) + cimag(a[1 * 3 + 0]));

    a_components[3] = (creal(a[0 * 3 + 0]) - creal(a[1 * 3 + 1]));

    a_components[4] = (creal(a[0 * 3 + 2]) + creal(a[2 * 3 + 0]));

    a_components[5] = (-cimag(a[0 * 3 + 2]) + cimag(a[2 * 3 + 0]));

    a_components[6] = (creal(a[1 * 3 + 2]) + creal(a[2 * 3 + 1]));

    a_components[7] = (-cimag(a[1 * 3 + 2]) + cimag(a[2 * 3 + 1]));

    a_components[8] = (creal(a[0 * 3 + 0]) + creal(a[1 * 3 + 1]) - 2.0 * creal(a[2 * 3 + 2])) / pow(3.0, 0.5);
}