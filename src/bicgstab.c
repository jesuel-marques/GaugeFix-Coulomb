#include <SU3_ops.h>
#include <fields.h>
#include <geometry.h>
#include <stdbool.h>
#include <stdlib.h>
#include <tgmath.h>

#define ELEM_VEC_INV(v, position, dirac, color) ((v) + (((((position.pos[T_INDX]) * lattice_param.n_SPC + (position.pos[X_INDX])) * lattice_param.n_SPC + (position.pos[Y_INDX])) * lattice_param.n_SPC + (position.pos[Z_INDX])) * 4 + (dirac)) * Nc + (color))

double complex Gamma[DIM][4][4];  // Matrizes gama de Dirac
// double complex Sigma[d][d][4][4];	//	Matrizes sigma de Dirac/* Vou ignorar improvements por enquanto */

double complex identityDirac[4][4];  // Identidade no espaço de Dirac

double complex identityDiracminusGamma[DIM][4][4];  //	Identidade menos gama de Dirac
double complex identityDiracplusGamma[DIM][4][4];   //	Identidade mais gama de Dirac

#define LOOP_DIRAC(alpha) for (alpha = 0; alpha < 4; alpha++)

static void copy_vector(Scalar *f, Scalar *copy_f, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        copy_f[i] = f[i];
    }
}

static void prodwithscalarVector(double complex num, Scalar *f, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = num * f[i];
    }
}

static void accumulateVector(Scalar *f, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] += f[i];
    }
}

static void subtractVector(Scalar *f, Scalar *g, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = f[i] - g[i];
    }
}

static void fmaVector(Scalar *f, Scalar num, Scalar *g, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = fma(num, g[i], f[i]);
    }
}

static void fmaaccumulateVector(Scalar *f, Scalar num, Scalar *g, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        f[i] += fma(num, g[i], f[i]);
    }
}

static void doublefmaVector(Scalar *f, Scalar num1, Scalar *g1, Scalar num2, Scalar *g2, Scalar *result, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        result[i] = fma(num2, g2[i], fma(num1, g1[i], f[i]));
    }
}

static void doublefmaaccumulateVector(Scalar *f, Scalar num1, Scalar *g1, Scalar num2, Scalar *g2, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        f[i] = fma(num2, g2[i], fma(num1, g1[i], f[i]));
    }
}

static Scalar dotproductVector(Scalar *f, Scalar *g, size_t number_of_elements) {
    Scalar dot_product;
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
                                prod1 = -kappa * identityDiracminusGamma[mu][alpha][alphap] * (link1.m[ELEM_3X3(a, ap)]) * (*(ELEM_VEC_INV(f, positionp1, alphap, ap)));

                                prod1 = -kappa * identityDiracplusGamma[mu][alpha][alphap] * (link2.m[ELEM_3X3(a, ap)]) * (*(ELEM_VEC_INV(f, positionp2, alphap, ap)));

                                //	Anti-periodic boundary condition
                                (*(ELEM_VEC_INV(g, position, alpha, a))) += (mu != T_INDX || position.pos[T_INDX] != (lattice_param.n_T - 1)) ? prod1 : -prod1;
                                (*(ELEM_VEC_INV(g, position, alpha, a))) += (mu != T_INDX || position.pos[T_INDX] != 0) ? prod2 : -prod2;
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
                    *(ELEM_VEC_INV(source, position, alpha, a)) = (samePositionQ(position, source_position) || alpha == dirac_index || a == color_index) ? 0.0 : 1.0;
                }
            }
        }
    }
}

void BiCGStab(int dirac_index, int color_index, double kappa, double tolerance, Mtrx3x3 *U, Scalar *source, Scalar *inverse_column) {
    //	Algoritmo para inversão da matriz de Dirac (descrito em Gattringer pg. 140 e Leandro pg. 110)

    // Variáveis do algoritmo de inversão

    Scalar *x;

    Scalar *r;
    Scalar *rtilde;

    Scalar *v;
    Scalar *p;
    Scalar *aux1;  // b e t na notação do Leandro
    Scalar *aux2;  // s na notação do Leandro

    const size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    x = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    r = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    rtilde = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    v = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    p = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    aux1 = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));
    aux2 = (Scalar *)calloc(sizeof_vector, sizeof(Scalar));

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

    copy_vector(source, aux1, sizeof_vector);

    copy_vector(aux1, x, sizeof_vector); /*x0=b*/

    actDiracOperator(U, kappa, x, aux2);          /*aux2=A.x0*/
    subtractVector(aux1, aux2, r, sizeof_vector); /*r=b-A.x0*/

    copy_vector(r, rtilde, sizeof_vector); /*rtilde=r*/

    do {
        rho = dotproductVector(rtilde, r, sizeof_vector); /*rho=rtildedag.r*/

        if (rho == 0.0) {
            printf("Method failed.");
        }

        if (first_iteration) {
            copy_vector(r, p, sizeof_vector); /*p=r*/
            first_iteration = false;
        }

        else {
            beta = alpha * rho / (omega * previous_rho);
            doublefmaVector(r, beta, p, -beta * omega, v, aux1, sizeof_vector); /*aux1=r+beta*(p-omega*v)*/
            copy_vector(aux1, p, sizeof_vector);                                /*p=r+beta*(p-omega*v)*/
        }
        actDiracOperator(U, kappa, p, v); /*v=A.p*/
        alpha = rho / dotproductVector(rtilde, v, sizeof_vector);
        fmaVector(r, -alpha, v, aux2, sizeof_vector); /*s=r-alfa*v*/

        /*auxa=sdag.s*/
        norm = fabs(creal(dotproductVector(aux2, aux2, sizeof_vector)));

        if (cont % 5 == 0)
            printf("norm s: %3.2E in Dirac %d Color %d\n", norm, dirac_index, color_index);

        if (norm < tolerance) {
            fmaaccumulateVector(x, alpha, p, sizeof_vector); /*x=x+alfa*p*/

            inverted = true;
        }

        else {
            actDiracOperator(U, kappa, aux2, aux1); /*t=A.s*/

            omega = dotproductVector(aux1, aux2, sizeof_vector) / dotproductVector(aux1, aux1, sizeof_vector); /*omega=tdag.s/tdag.t*/

            fmaVector(aux2, -omega, aux1, r, sizeof_vector);

            doublefmaaccumulateVector(x, alpha, p, omega, aux2, sizeof_vector);
            /*x=x+alfa*p+omega*s*/

            previous_rho = rho;
        }

        cont++;

    } while (!inverted);

    CopiarVetorInversao(x, inverse_column);

    free(x);

    free(r);
    free(rtilde);

    free(v);
    free(p);

    free(aux1);
    free(aux2);
}
