#include <SU3_ops.h>
#include <dirac.h>
#include <fields.h>
#include <geometry.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#define ELEM_VEC_INV(v, position, dirac, color) ((v) + (((((position.pos[T_INDX]) * lattice_param.n_SPC + (position.pos[X_INDX])) * lattice_param.n_SPC + (position.pos[Y_INDX])) * lattice_param.n_SPC + (position.pos[Z_INDX])) * 4 + (dirac)) * Nc + (color))

DiracMatrix gamma[DIM] = {{{{0, 0, 1, 0},
                            {0, 0, 0, 1},
                            {1, 0, 0, 0},
                            {0, 1, 0, 0}}},

                          {{{0, 0, 0, -I},
                            {0, 0, -I, 0},
                            {0, I, 0, 0},
                            {I, 0, 0, 0}}},

                          {{{0, 0, 0, -1},
                            {0, 0, 1, 0},
                            {0, 1, 0, 0},
                            {-1, 0, 0, 0}}},

                          {{{0, 0, -I, 0},
                            {0, 0, 0, I},
                            {I, 0, 0, 0},
                            {0, -I, 0, 0}}}};

DiracMatrix identity_plus_gamma[DIM] = {{{{1, 0, 1, 0},
                                          {0, 1, 0, 1},
                                          {1, 0, 1, 0},
                                          {0, 1, 0, 1}}},

                                        {{{1, 0, 0, -I},
                                          {0, 1, -I, 0},
                                          {0, I, 1, 0},
                                          {I, 0, 0, 1}}},

                                        {{{1, 0, 0, -1},
                                          {0, 1, 1, 0},
                                          {0, 1, 1, 0},
                                          {-1, 0, 0, 1}}},

                                        {{{1, 0, -I, 0},
                                          {0, 1, 0, I},
                                          {I, 0, 1, 0},
                                          {0, -I, 0, 1}}}};

static DiracMatrix identity_minus_gamma[4] = {{{{1, 0, -1, 0},
                                                {0, 1, 0, -1},
                                                {-1, 0, 1, 0},
                                                {0, -1, 0, 1}}},

                                              {{{1, 0, 0, I},
                                                {0, 1, I, 0},
                                                {0, -I, 1, 0},
                                                {-I, 0, 0, 1}}},

                                              {{{1, 0, 0, 1},
                                                {0, 1, -1, 0},
                                                {0, -1, 1, 0},
                                                {1, 0, 0, 1}}},

                                              {{{1, 0, I, 0},
                                                {0, 1, 0, -I},
                                                {-I, 0, 1, 0},
                                                {0, I, 0, 1}}}};

static DiracMatrix identityDirac = {{{1, 0, 0, 0},
                                     {0, 1, 0, 0},
                                     {0, 0, 1, 0},
                                     {0, 0, 0, 1}}};

static void copy_vector(Scalar *f, Scalar *copy_f, size_t number_of_elements) {
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

void actDiracOperator(Mtrx3x3 *U, double kappa, Scalar *f, Scalar *g) {
    //	Calcula g=D.f, ou seja, a matriz de Dirac agindo no vetor f, e coloca o resultado em g.

    PosVec position;
    PosVec positionp1;
    PosVec positionp2;

    Mtrx3x3 link1;
    Scalar prod1;

    Mtrx3x3 link2;
    Scalar prod2;

    LorentzIdx t;
    MtrxIdx3 a, ap;
    int alpha, alphap;
    LorentzIdx mu;

    LOOP_TEMPORAL(t) {
        position.pos[T_INDX] = t;
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    // Diagonal term: simply copy content of g to f

                    (*(ELEM_VEC_INV(g, position, alpha, a))) = (*(ELEM_VEC_INV(f, position, alpha, a)));

                    // Hopping term with forward (positive mu, 1) and backward
                    // (negative mu, 2) neighbors

                    LOOP_LORENTZ(mu) {
                        getLinkMatrix(U, position, mu, FRONT, &link1);
                        getLinkMatrix(U, position, mu, REAR, &link2);

                        positionp1 = getNeighbour(position, mu, FRONT);
                        positionp2 = getNeighbour(position, mu, REAR);

                        LOOP_DIRAC(alphap) {
                            LOOP_3(ap) {
                                prod1 = -kappa *
                                        identity_minus_gamma[mu].m[alpha][alphap] *
                                        (link1.m[ELEM_3X3(a, ap)]) *
                                        (*(ELEM_VEC_INV(f, positionp1, alphap, ap)));

                                prod2 = -kappa *
                                        identity_plus_gamma[mu].m[alpha][alphap] *
                                        (link2.m[ELEM_3X3(a, ap)]) *
                                        (*(ELEM_VEC_INV(f, positionp2, alphap, ap)));

                                //	Anti-periodic boundary condition
                                (*(ELEM_VEC_INV(g, position, alpha, a))) +=
                                    (mu != T_INDX ||
                                     position.pos[T_INDX] != (lattice_param.n_T - 1))
                                        ? prod1
                                        : -prod1;
                                (*(ELEM_VEC_INV(g, position, alpha, a))) +=
                                    (mu != T_INDX ||
                                     position.pos[T_INDX] != 0)
                                        ? prod2
                                        : -prod2;
                            }
                        }
                    }
                }
            }
        }
    }
    // Melhoria de trevo /* Vou ignorar improvements por enquanto */

    // for(int alfal=0;alfal<4;alfal++)
    // 	for(int al=0;al<3;al++){
    // 		(*(ELEM_VEC_INV(g, posicao, alfa, a))) += cSWkappaSigmaF[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][alfal][al] * (*(ELEM_VEC_INV(f, posicao, alfal, al)));
    // }
}

void initializePointSource(PosVec source_position, int dirac_index, int color_index, double complex *source) {
    //	Inicializa b, que é um vetor coluna da matriz obtida pelas deltas de posição, índices de Dirac e cor

    LorentzIdx t;
    PosVec position;
    MtrxIdx3 a;
    int alpha;
    LOOP_TEMPORAL(t) {
        position.pos[T_INDX] = t;
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    *(ELEM_VEC_INV(source, position, alpha, a)) =
                        (samePositionQ(position, source_position) &&
                         alpha == dirac_index &&
                         a == color_index)
                            ? 1.0
                            : 0.0;
                    // if (*(ELEM_VEC_INV(source, position, alpha, a)) == 1.0) {
                    //     printf("dirac: %d color: %d \n", dirac_index, color_index);
                    //     printPosVec(position);
                    //     getchar();
                }
            }
        }
    }
}

void BiCGStab(double kappa, Mtrx3x3 *U, Scalar *source, Scalar *inverse_column, double tolerance) {
    //	Algoritmo para inversão da matriz de Dirac (descrito em Gattringer pg. 140 e Leandro pg. 110)

    // Variáveis do algoritmo de inversão

    Scalar *x;

    Scalar *r;
    Scalar *rtilde;

    Scalar *v;
    Scalar *p;
    Scalar *t;  // b e t na notação do Leandro
    Scalar *s;  // s na notação do Leandro

    const size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    x = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    r = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    rtilde = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    v = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    p = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    t = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    s = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

    int cont = 0;

    Scalar auxa;
    Scalar auxb;

    Scalar alpha;
    Scalar beta;
    Scalar omega;

    Scalar rho, previous_rho;

    double norm;

    bool first_iteration = true;
    bool inverted = false;

    // copy_vector(source, t, sizeof_vector);

    copy_vector(source, x, sizeof_vector); /*x0=b*/

    actDiracOperator(U, kappa, x, s); /*s=A.x0*/
    // printf("%.18lf\n", *ELEM_VEC_INV(s, assignPosition(2, 1, 4, 3), 1, 2));
    subtractVector(source, s, r, sizeof_vector); /*r=b-A.x0*/

    copy_vector(r, rtilde, sizeof_vector); /*rtilde=r*/

    do {
        rho = dotproductVector(rtilde, r, sizeof_vector); /*rho=rtildedag.r*/
        // printf("alpha: %lf beta: %lf omega: %lf rho: %lf\n", alpha, beta, omega, rho);
        // getchar();
        if (rho == 0.0) {
            printf("Method failed.");
        }

        if (first_iteration) {
            copy_vector(r, p, sizeof_vector); /*p=r*/
            first_iteration = false;
        }

        else {
            beta = alpha * rho / (omega * previous_rho);
            doublefmaVector(r, beta, p, -beta * omega, v, t, sizeof_vector); /*t=r+beta*(p-omega*v)*/
            copy_vector(t, p, sizeof_vector);                                /*p=r+beta*(p-omega*v)*/
        }
        actDiracOperator(U, kappa, p, v); /*v=A.p*/
        alpha = rho / dotproductVector(rtilde, v, sizeof_vector);
        fmaVector(r, -alpha, v, s, sizeof_vector); /*s=r-alpha*v*/

        /*auxa=sdag.s*/
        norm = fabs(creal(dotproductVector(s, s, sizeof_vector)));

        if (cont % 5 == 0)
            printf("norm s: %3.2E\n", norm);

        if (norm < tolerance) {
            fmaaccumulateVector(x, alpha, p, sizeof_vector); /*x=x+alfa*p*/

            inverted = true;
        }

        else {
            actDiracOperator(U, kappa, s, t); /*t=A.s*/

            omega = dotproductVector(t, s, sizeof_vector) / dotproductVector(t, t, sizeof_vector); /*omega=tdag.s/tdag.t*/

            fmaVector(s, -omega, t, r, sizeof_vector); /* r = s - omega*t*/

            doublefmaaccumulateVector(x, alpha, p, omega, s, sizeof_vector);
            /*x=x+alfa*p+omega*s*/

            previous_rho = rho;
        }

        cont++;

    } while (!inverted);

    copy_vector(x, inverse_column, sizeof_vector);

    free(x);

    free(r);
    free(rtilde);

    free(v);
    free(p);

    free(t);
    free(s);
}
