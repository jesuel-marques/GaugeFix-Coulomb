#include <bicgstab.h>
#include <stdbool.h>
#include <stdlib.h>

static void copyVector(Scalar *f, Scalar *copy_f, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        copy_f[i] = f[i];
    }
}

static void prodwithscalarVector(Scalar num, Scalar *f, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = num * f[i];
    }
}

static void subtractVector(Scalar *f, Scalar *g, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = f[i] - g[i];
    }
}

static void fmaVector(Scalar *f, Scalar num, Scalar *g, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = f[i] + num * g[i];
    }
}

static void fmaaccumulateVector(Scalar *f, Scalar num, Scalar *g, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        f[i] += num * g[i];
    }
}

static void doublefmaVector(Scalar *f, Scalar num1, Scalar *g1, Scalar num2, Scalar *g2, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = f[i] + num1 * g1[i] + num2 * g2[i];
    }
}

static void doublefmaaccumulateVector(Scalar *f, Scalar num1, Scalar *g1, Scalar num2, Scalar *g2, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        f[i] += num1 * g1[i] + num2 * g2[i];
    }
}

static Scalar dotproductVector(Scalar *f, Scalar *g, size_t number_of_elements) {
    Scalar dot_product = 0.0;
    for (size_t i = 0; i < number_of_elements; i++) {
        dot_product += conj(f[i]) * g[i];
    }
    return dot_product;
}

double BiCGStab(void (*operator)(Scalar *, Scalar *), Scalar *source, Scalar *inverse_column, const double tolerance, const size_t sizeof_vector) {
    //	Algoritmo para inversão da matriz de Dirac (descrito em Gattringer pg. 140 e Leandro pg. 110)

    // Variáveis do algoritmo de inversão

    Scalar *x = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar *r = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *rtilde = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar *v = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *p = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *t = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));  // b e t na notação do Leandro
    Scalar *s = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));  // s na notação do Leandro

    Scalar auxa, auxb;
    Scalar alpha, beta, omega;
    Scalar rho, previous_rho;

    double norm;

    bool first_iteration = true;
    bool inverted = false;

    // copyVector(source, t, sizeof_vector);

    copyVector(source, x, sizeof_vector); /*x0=b*/

    operator(x, s);                              /*s=A.x0*/
    subtractVector(source, s, r, sizeof_vector); /*r=b-A.x0*/

    copyVector(r, rtilde, sizeof_vector); /*rtilde=r*/

    do {
        rho = dotproductVector(rtilde, r, sizeof_vector); /*rho=rtildedag.r*/
        // printf("alpha: %lf beta: %lf omega: %lf rho: %lf\n", alpha, beta, omega, rho);
        // getchar();
        if (rho == 0.0) {
            printf("Method failed.");
        }

        if (first_iteration) {
            copyVector(r, p, sizeof_vector); /*p=r*/
            first_iteration = false;
        }

        else {
            beta = alpha * rho / (omega * previous_rho);
            doublefmaVector(r, beta, p, -beta * omega, v, t, sizeof_vector); /*t=r+beta*(p-omega*v)*/
            copyVector(t, p, sizeof_vector);                                 /*p=r+beta*(p-omega*v)*/
        }
        operator(p, v); /*v=A.p*/
        alpha = rho / dotproductVector(rtilde, v, sizeof_vector);
        fmaVector(r, -alpha, v, s, sizeof_vector); /*s=r-alpha*v*/

        /*auxa=sdag.s*/
        norm = creal(dotproductVector(s, s, sizeof_vector));

        // if (cont % 5 == 0)
        //     printf("norm s: %3.2E\n", norm);

        if (norm < tolerance) {
            fmaaccumulateVector(x, alpha, p, sizeof_vector); /*x=x+alfa*p*/

            inverted = true;
        }

        else {
            operator(s, t); /*t=A.s*/

            omega = dotproductVector(t, s, sizeof_vector) / dotproductVector(t, t, sizeof_vector); /*omega=tdag.s/tdag.t*/

            fmaVector(s, -omega, t, r, sizeof_vector); /* r = s - omega*t*/

            doublefmaaccumulateVector(x, alpha, p, omega, s, sizeof_vector);
            /*x=x+alfa*p+omega*s*/

            previous_rho = rho;
        }

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
