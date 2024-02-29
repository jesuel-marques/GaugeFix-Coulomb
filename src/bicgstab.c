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

    return f;
}

static Scalar *subtractVector(Scalar *restrict f, Scalar *restrict g, Scalar *restrict result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = f[i] - g[i];
    }
    return result;
}

static Scalar *fmaVector(Scalar *f, Scalar num, Scalar *g, Scalar *result, size_t number_of_elements) {
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

static Scalar *doublefmaVector(Scalar *f, Scalar num1, Scalar *g1, Scalar num2, Scalar *g2, Scalar *result, size_t number_of_elements) {
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

static Scalar *productwithScalarVector(Scalar num, Scalar *f, Scalar *numf, size_t number_of_elements) {
    Scalar dot_product = 0.0;

    for (size_t i = 0; i < number_of_elements; i++) {
        numf[i] = num * f[i];
    }
    return numf;
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

double BiCGStabM(Scalar *(*operator)(Scalar *, Scalar *), Scalar *restrict source, Scalar *restrict inverse_column, const double sigma, const double tolerance, const size_t sizeof_vector) {
    //	Multimass version of BiCGStab by Beat Jegerlehner (arXiv:hep-lat/9612014v1)
    bool first_iteration = true;

    Scalar *x_sigma = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar *old_r = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *r = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *s = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *As = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar *s_sigma = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar *w = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *w0 = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *Aw = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar alpha = 0.0, old_beta = 1.0, beta, old_delta, delta, phi, chi;
    Scalar alpha_sigma = 0.0, beta_sigma, old_zeta_sigma = 1.0, twice_old_zeta_sigma = 0, zeta_sigma, old_rho_sigma = 1.0, rho_sigma, chi_sigma;

    double norm_w;

    setVector(x_sigma, 0.0, sizeof_vector);
    copyVector(source, old_r, sizeof_vector);
    copyVector(source, s, sizeof_vector);
    copyVector(source, s_sigma, sizeof_vector);

    setRandomVector(w0, sizeof_vector);
    copyVector(w0, w, sizeof_vector);
    operator(s, As);
    while (true) {
        if (first_iteration == true) {
            first_iteration = false;
            old_delta = dotproductVector(w0, r, sizeof_vector);

            phi = dotproductVector(w, As, sizeof_vector) / delta;
        }
        beta = -1 / phi;

        zeta_sigma = (old_zeta_sigma * twice_old_zeta_sigma * old_beta) / (beta * alpha * (twice_old_zeta_sigma - old_zeta_sigma) + twice_old_zeta_sigma * old_beta * (1 - sigma * beta));

        beta_sigma = beta * zeta_sigma / old_zeta_sigma;

        fmaVector(old_r, beta, As, w, sizeof_vector);
        if ((norm_w = dotproductVector(w, w, sizeof_vector)) < tolerance) {
            // I'm guessing here because the paper doesn't explain the stop criterion
            fmaaccumulateVector(x_sigma, -beta_sigma, s, sizeof_vector);
            break;
        }
        operator(w, Aw);

        chi = dotproductVector(Aw, w, sizeof_vector) / dotproductVector(Aw, Aw, sizeof_vector);
        chi_sigma = chi / (1 + chi * sigma);
        rho_sigma = old_rho_sigma / (1 + chi * sigma);

        fmaVector(w, -chi, Aw, r, sizeof_vector);

        doublefmaaccumulateVector(x_sigma, -beta_sigma, s, chi_sigma * rho_sigma * zeta_sigma, w, sizeof_vector);

        delta = dotproductVector(w0, r, sizeof_vector);

        alpha = -beta * delta / (old_delta * chi);
        alpha_sigma = alpha * zeta_sigma * beta_sigma / (old_zeta_sigma * beta);

        doublefmaVector(r, alpha, s, -alpha * chi, As, s, sizeof_vector);
        operator(s, As);

        productwithScalarVector(zeta_sigma * rho_sigma, r, s_sigma, sizeof_vector);
        fmaaccumulateVector(s_sigma, alpha_sigma, s_sigma, sizeof_vector);
        fmaaccumulateVector(s_sigma, -alpha_sigma * chi_sigma * zeta_sigma * rho_sigma / beta_sigma, w, sizeof_vector);
        fmaaccumulateVector(s_sigma, alpha_sigma * (chi_sigma / beta_sigma) * old_zeta_sigma * old_rho_sigma, old_r, sizeof_vector);

        phi = dotproductVector(w0, As, sizeof_vector) / delta;

        old_beta = beta;
        old_delta = delta;
        twice_old_zeta_sigma = old_zeta_sigma;
        old_zeta_sigma = zeta_sigma;
        old_rho_sigma = rho_sigma;

        copyVector(r, old_r, sizeof_vector);
    }
    return norm_w;
}

double CG(Scalar *(*operator)(Scalar *, Scalar *), Scalar *restrict source, Scalar *restrict inverse_column, const double tolerance, const size_t sizeof_vector) {
    //	Algoritmo para inversão da matriz (descrito em mscg.pdf do Martin Lüscher)

    // Variáveis do algoritmo de inversão

    Scalar *x = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *r = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *p = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    Scalar *Ap = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    Scalar alpha, beta;

    Scalar old_norm_r = 0.0, norm_r = 0.0;

    setVector(x, 0.0, sizeof_vector);
    copyVector(source, r, sizeof_vector);
    copyVector(source, p, sizeof_vector);

    old_norm_r = dotproductVector(r, r, sizeof_vector);
    printf("norm r: %3.2E\n", creal(old_norm_r));
    system("sleep 5");
    int cont = 0;
    do {
        operator(p, Ap);
        alpha = old_norm_r / dotproductVector(p, Ap, sizeof_vector);
        fmaaccumulateVector(x, alpha, p, sizeof_vector);

        fmaaccumulateVector(r, -alpha, Ap, sizeof_vector);

        norm_r = dotproductVector(r, r, sizeof_vector);

        printf("norm r: %3.2E\n", creal(norm_r));

        if (creal(norm_r) < tolerance) {
            break;
        } else {
            beta = norm_r / old_norm_r;
            fmaVector(r, beta, p, p, sizeof_vector);
            old_norm_r = norm_r;
        }

        cont++;
    } while (true);

    copyVector(x, inverse_column, sizeof_vector);

    free(x);

    free(r);

    free(p);
    free(Ap);

    return norm_r;
}
