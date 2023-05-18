#include <SU3_ops.h>
#include <bicgstab.h>
#include <dirac.h>
#include <tgmath.h>

#define ELEM_VEC_INV(v, position, dirac, color) ((v) + (((((position.pos[T_INDX]) * lattice_param.n_SPC + (position.pos[Z_INDX])) * lattice_param.n_SPC + (position.pos[Y_INDX])) * lattice_param.n_SPC + (position.pos[X_INDX])) * 4 + (dirac)) * Nc + (color))

static Mtrx3x3 *U_dirac_operator;
static double kappa_dirac_operator;

static DiracMatrix gamma[DIM] = {
    {{{0, 0, 1, 0},
      {0, 0, 0, 1},
      {1, 0, 0, 0},
      {0, 1, 0, 0}}},

    {{{0, 0, -I, 0},
      {0, 0, 0, I},
      {I, 0, 0, 0},
      {0, -I, 0, 0}}},

    {{{0, 0, 0, -1},
      {0, 0, 1, 0},
      {0, 1, 0, 0},
      {-1, 0, 0, 0}}},

    {{{0, 0, 0, -I},
      {0, 0, -I, 0},
      {0, I, 0, 0},
      {I, 0, 0, 0}}}

};

static DiracMatrix identity_plus_gamma[DIM] = {{{{1, 0, 1, 0},
                                                 {0, 1, 0, 1},
                                                 {1, 0, 1, 0},
                                                 {0, 1, 0, 1}}},

                                               {{{1, 0, -I, 0},
                                                 {0, 1, 0, I},
                                                 {I, 0, 1, 0},
                                                 {0, -I, 0, 1}}},

                                               {{{1, 0, 0, -1},
                                                 {0, 1, 1, 0},
                                                 {0, 1, 1, 0},
                                                 {-1, 0, 0, 1}}},

                                               {{{1, 0, 0, -I},
                                                 {0, 1, -I, 0},
                                                 {0, I, 1, 0},
                                                 {I, 0, 0, 1}}}};

static DiracMatrix identity_minus_gamma[DIM] = {{{{1, 0, -1, 0},
                                                  {0, 1, 0, -1},
                                                  {-1, 0, 1, 0},
                                                  {0, -1, 0, 1}}},

                                                {{{1, 0, I, 0},
                                                  {0, 1, 0, -I},
                                                  {-I, 0, 1, 0},
                                                  {0, I, 0, 1}}},

                                                {{{1, 0, 0, 1},
                                                  {0, 1, -1, 0},
                                                  {0, -1, 1, 0},
                                                  {1, 0, 0, 1}}},

                                                {{{1, 0, 0, I},
                                                  {0, 1, I, 0},
                                                  {0, -I, 1, 0},
                                                  {-I, 0, 0, 1}}}};

static DiracMatrix identityDirac = {{{1, 0, 0, 0},
                                     {0, 1, 0, 0},
                                     {0, 0, 1, 0},
                                     {0, 0, 0, 1}}};

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

void DiracOperator(Scalar *f, Scalar *g) {
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

                    (*(ELEM_VEC_INV(g, position, alpha, a))) =
                        (*(ELEM_VEC_INV(f, position, alpha, a)));

                    // Hopping term with forward (positive mu, 1) and backward
                    // (negative mu, 2) neighbors

                    LOOP_LORENTZ(mu) {
                        getLinkMatrix(U_dirac_operator, position, mu, FRONT, &link1);
                        getLinkMatrix(U_dirac_operator, position, mu, REAR, &link2);

                        positionp1 = getNeighbour(position, mu, FRONT);
                        positionp2 = getNeighbour(position, mu, REAR);

                        LOOP_DIRAC(alphap) {
                            LOOP_3(ap) {
                                prod1 = -kappa_dirac_operator *
                                        identity_minus_gamma[mu].m[alpha][alphap] *
                                        (link1.m[ELEM_3X3(a, ap)]) *
                                        (*(ELEM_VEC_INV(f, positionp1, alphap, ap)));

                                prod2 = -kappa_dirac_operator *
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
                    // Melhoria de trevo /* Vou ignorar improvements por enquanto */

                    // for(int alfal=0;alfal<4;alfal++)
                    // 	for(int al=0;al<3;al++){
                    // 		(*(ELEM_VEC_INV(g, posicao, alfa, a))) += cSWkappaSigmaF[posicao[0]][posicao[1]][posicao[2]][posicao[3]][alfa][a][alfal][al] * (*(ELEM_VEC_INV(f, posicao, alfal, al)));
                    // }
                }
            }
        }
    }
}

double invertDiracOperator(double kappa, Mtrx3x3 *U, Scalar *source, Scalar *inverse_column, double tolerance, double (*inversion_algorithm)(void (*)(Scalar *, Scalar *), Scalar *, Scalar *, double, size_t)) {
    U_dirac_operator = U;
    kappa_dirac_operator = kappa;

    size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    return inversion_algorithm(DiracOperator, source, inverse_column, tolerance, sizeof_vector);
}

void printDiracOP(Mtrx3x3 *U, double kappa, FILE *file_dirac_op) {
    //	Imprime a Matriz de Dirac explicitamente (Comparação com Mathematica)

    PosVec position1, position2;

    Scalar *f;
    Scalar *g;

    f = (double complex *)calloc(lattice_param.volume * 4 * Nc, sizeof(double complex));
    g = (double complex *)calloc(lattice_param.volume * 4 * Nc, sizeof(double complex));

    int alpha, beta;
    MtrxIdx3 a, b;
    LorentzIdx t1, t2;
    LOOP_TEMPORAL(t1) {
        position1.pos[T_INDX] = t1;
        LOOP_SPATIAL(position1) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    (*(ELEM_VEC_INV(f, position1, alpha, a))) = 1.0;
                    DiracOperator(f, g);
                    (*(ELEM_VEC_INV(f, position1, alpha, a))) = 0.0;
                    LOOP_TEMPORAL(t2) {
                        position2.pos[T_INDX] = t2;
                        LOOP_SPATIAL(position2) {
                            LOOP_DIRAC(beta) {
                                LOOP_3(b) {
                                    fprintf(file_dirac_op, "%.20lf+I(%.20lf)\t",
                                            creal((*(ELEM_VEC_INV(g,
                                                                  position2,
                                                                  beta,
                                                                  b)))),
                                            cimag((*(ELEM_VEC_INV(g,
                                                                  position2,
                                                                  beta,
                                                                  b)))));
                                }
                            }
                        }
                    }
                    fprintf(file_dirac_op, "\n");
                }
            }
        }
    }

    free(f);
    free(g);
}

void printInverse(char *inverse_filename, Scalar *inverse[4][Nc]) {
    //	Imprime em um arquivo as colunas da matriz inversa do operador de Dirac

    PosVec position;

    size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    FILE *inverse_file;
    if ((inverse_file = fopen(inverse_filename, "wb")) == NULL) {
        fprintf(stderr,
                "Error: Problem creating file %s for inverse of Dirac matrix.\n", inverse_filename);
    }

    int alpha;
    MtrxIdx3 a;
    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            if (fwrite(inverse[alpha][a],
                       sizeof(Scalar),
                       sizeof_vector,
                       inverse_file) != sizeof_vector) {
                fprintf(stderr, "Error: Problem writing to file %s for inverse of Dirac matrix.\n", inverse_filename);
            }
        }
    }
    fclose(inverse_file);
}