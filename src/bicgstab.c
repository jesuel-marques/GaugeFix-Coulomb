#include <bicgstab.h>
#include <ranlux.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static Scalar *copyVector(Scalar *restrict f, Scalar *restrict copy_f, size_t number_of_elements) {
    memcpy(copy_f, f, number_of_elements * sizeof(Scalar));
    return copy_f;
}

static Scalar *setVector(Scalar *restrict f, double complex value, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        f[i] = value;
    }
    return f;
}

static Scalar *setRandomVector(Scalar *restrict f, size_t number_of_elements) {
    rlxd_init(2, 32244000);
    ranlxd(f, 2 * number_of_elements);
}

static Scalar *subtractVector(Scalar *restrict f, Scalar *restrict g, Scalar *restrict result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = f[i] - g[i];
    }
    return result;
}

static Scalar *fmaVector(Scalar *restrict f, Scalar num, Scalar *restrict g, Scalar *restrict result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = f[i] + num * g[i];
    }
    return result;
}

static Scalar *fmaaccumulateVector(Scalar *restrict f, Scalar num, Scalar *restrict g, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        f[i] += num * g[i];
    }
    return f;
}

static Scalar *doublefmaVector(Scalar *restrict f, Scalar num1, Scalar *restrict g1, Scalar num2, Scalar *restrict g2, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = f[i] + num1 * g1[i] + num2 * g2[i];
    }

    return result;
}

static Scalar *doublefmaaccumulateVector(Scalar *restrict f, Scalar num1, Scalar *restrict g1, Scalar num2, Scalar *restrict g2, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        f[i] += num1 * g1[i] + num2 * g2[i];
    }
    return f;
}

static Scalar dotproductVector(Scalar *f, Scalar *g, size_t number_of_elements) {
    Scalar dot_product = 0.0;

    for (size_t i = 0; i < number_of_elements; i++) {
        dot_product += conj(f[i]) * g[i];
    }
    return dot_product;
}

double BiCGStab(Scalar *(*operator)(Scalar *, Scalar *), Scalar *restrict source, Scalar *restrict inverse_column, const double tolerance, const size_t sizeof_vector) {
    //	Algoritmo para inversão da matriz de Dirac (descrito em Gattringer pg. 140 e Leandro pg. 110)

    // Variáveis do algoritmo de inversão

    Scalar *x = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar *r = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *rtilde = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar *v = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *p = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *t = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *s = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar alpha, beta, omega;
    Scalar rho, previous_rho;

    double norm;

    bool first_iteration = true;
    bool inverted = false;

    // copyVector(source, t, sizeof_vector);

    // copyVector(source, x, sizeof_vector); /*x0=b*/

    setRandomVector(x, sizeof_vector);
    // setVector(x, 1.0, sizeof_vector);

    /*s=A.x0*/
    // subtractVector(source, s, r, sizeof_vector); /*r=b-A.x0*/

    copyVector(subtractVector(source, operator(x, s), r, sizeof_vector), rtilde, sizeof_vector); /*rtilde=r*/
    // copyVector(x, rtilde, sizeof_vector); /* rtilde = x*/

    int cont = 0;
    do {
        rho = dotproductVector(rtilde, r, sizeof_vector); /*rho=rtildedag.r*/
        if (rho == 0.0) {
            printf("Method failed.");
        }

        if (first_iteration) {
            copyVector(r, p, sizeof_vector); /*p=r*/
            first_iteration = false;
        }

        else {
            beta = alpha * rho / (omega * previous_rho);
            doublefmaVector(r, beta, p, -beta * omega, v, t, sizeof_vector);
            /*t=r+beta*(p-omega*v) => t = r + beta * p - beta * omega * v */
            copyVector(t, p, sizeof_vector);
            /*p=r+beta*(p-omega*v) => p = r + beta * p - beta * omega * v */
        }

        /*v=A.p*/
        alpha = rho / dotproductVector(rtilde, operator(p, v), sizeof_vector);
        fmaVector(r, -alpha, v, s, sizeof_vector); /*s=r-alpha*v*/

        /*auxa=sdag.s*/
        norm = creal(dotproductVector(s, s, sizeof_vector));
        printf("norm s: %3.2E\n", norm);

        if (norm < tolerance) {
            fmaaccumulateVector(x, alpha, p, sizeof_vector); /*x=x+alfa*p*/

            inverted = true;
        } else {
            /*t=A.s*/

            operator(s, t);
            omega = dotproductVector(t, s, sizeof_vector) / dotproductVector(t, t, sizeof_vector); /*omega=tdag.s/tdag.t*/

            fmaVector(s, -omega, t, r, sizeof_vector); /* r = s - omega*t*/

            doublefmaaccumulateVector(x, alpha, p, omega, s, sizeof_vector);
            /*x+=alfa*p+omega*s*/

            previous_rho = rho;
        }

        cont++;
    } while (!inverted);

    copyVector(x, inverse_column, sizeof_vector);

    free(x);

    free(r);
    free(rtilde);

    free(v);
    free(p);

    free(t);
    free(s);

    return norm;
}
