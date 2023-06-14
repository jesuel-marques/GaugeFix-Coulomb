#include <SU3_ops.h>
#include <bicgstab.h>
#include <dirac.h>
#include <fields.h>
#include <plaquette.h>
#include <stdlib.h>
#include <tgmath.h>

static Mtrx3x3 *U_dirac_operator;
static DiracColorMatrix *sigmamunuFmunu_dirac_operator;

static double kappa_dirac_operator;
static double c_SW_dirac_operator;

static const DiracMatrix gamma[DIM] = {{{0, 0, 1, 0,
                                         0, 0, 0, 1,
                                         1, 0, 0, 0,
                                         0, 1, 0, 0}},

                                       {{0, 0, -I, 0,
                                         0, 0, 0, I,
                                         I, 0, 0, 0,
                                         0, -I, 0, 0}},

                                       {{0, 0, 0, -1,
                                         0, 0, 1, 0,
                                         0, 1, 0, 0,
                                         -1, 0, 0, 0}},

                                       {{0, 0, 0, -I,
                                         0, 0, -I, 0,
                                         0, I, 0, 0,
                                         I, 0, 0, 0}}};

// static DiracMatrix identity_plus_gamma[DIM] = {{{1, 0, 1, 0,
//                                                  0, 1, 0, 1,
//                                                  1, 0, 1, 0,
//                                                  0, 1, 0, 1}},  //    gamma[T_INDX]

//                                                {{1, 0, -I, 0,
//                                                  0, 1, 0, I,
//                                                  I, 0, 1, 0,
//                                                  0, -I, 0, 1}},  //    gamma[Z_INDX]

//                                                {{1, 0, 0, -1,
//                                                  0, 1, 1, 0,
//                                                  0, 1, 1, 0,
//                                                  -1, 0, 0, 1}},  //    gamma[Y_INDX]

//                                                {{1, 0, 0, -I,
//                                                  0, 1, -I, 0,
//                                                  0, I, 1, 0,
//                                                  I, 0, 0, 1}}};  //    gamma[X_INDX]

// static DiracMatrix identity_minus_gamma[DIM] = {{{1, 0, -1, 0,
//                                                   0, 1, 0, -1,
//                                                   -1, 0, 1, 0,
//                                                   0, -1, 0, 1}},

//                                                 {{1, 0, I, 0,
//                                                   0, 1, 0, -I,
//                                                   -I, 0, 1, 0,
//                                                   0, I, 0, 1}},

//                                                 {{1, 0, 0, 1,
//                                                   0, 1, -1, 0,
//                                                   0, -1, 1, 0,
//                                                   1, 0, 0, 1}},

//                                                 {{1, 0, 0, I,
//                                                   0, 1, I, 0,
//                                                   0, -I, 1, 0,
//                                                   -I, 0, 0, 1}}};

static const DiracMatrix identityDirac = {{1, 0, 0, 0,
                                           0, 1, 0, 0,
                                           0, 0, 1, 0,
                                           0, 0, 0, 1}};

#define DIGITS_MATRIX 8
void printDiracMatrix(DiracMatrix u) {
    unsigned short a, b;
    printf("{");
    LOOP_DIRAC(a) {
        printf("{");
        LOOP_DIRAC(b) {
            printf("%.*lf + I*(%.*lf)", DIGITS_MATRIX, creal(u.m[a * 4 + b]),
                   DIGITS_MATRIX, cimag(u.m[a * 4 + b]));

            b != 4 - 1 ? printf(", ") : 0;
        }
        a != 4 - 1 ? printf("},\n") : 0;
    }
    printf("}}\n\n");
}

DiracMatrix prodDirac(DiracMatrix m1, DiracMatrix m2) {
    DiracMatrix prod;
    unsigned short alpha, beta, gamma;
    LOOP_DIRAC(alpha) {
        LOOP_DIRAC(beta) {
            prod.m[alpha * 4 + beta] = 0.0;
            LOOP_DIRAC(gamma) {
                prod.m[alpha * 4 + beta] += m1.m[alpha * 4 + gamma] *
                                            m2.m[gamma * 4 + beta];
            }
        }
    }

    return prod;
}

void initializePointSource(PosVec source_position, DiracIdx dirac_index, MtrxIdx3 color_index, Scalar *source) {
    //	Inicializa b, que é um vetor coluna da matriz obtida pelas deltas de posição, índices de Dirac e cor

    PosVec position;
    MtrxIdx3 a;
    DiracIdx alpha;
    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    *(ELEM_VEC_POSDC(source, position, alpha, a)) =
                        (samePositionQ(position, source_position) &&
                         alpha == dirac_index &&
                         a == color_index)
                            ? 1.0
                            : 0.0;
                    // if (*(ELEM_VEC_POSDC(source, position, alpha, a)) == 1.0) {
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

    MtrxIdx3 a, ap;
    DiracIdx alpha, alphap;
    LorentzIdx mu;

    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    // Diagonal term: simply copy content of g to f

                    (*(ELEM_VEC_POSDC(g, position, alpha, a))) = (1.0 / (2.0 * kappa_dirac_operator)) * (*(ELEM_VEC_POSDC(f, position, alpha, a)));

                    // Hopping term with forward (positive mu, 1) and backward
                    // (negative mu, 2) neighbors

                    LOOP_LORENTZ(mu) {
                        getLinkMatrix(U_dirac_operator, position, mu, FRONT, &link1);
                        getLinkMatrix(U_dirac_operator, position, mu, REAR, &link2);

                        positionp1 = getNeighbour(position, mu, FRONT);
                        positionp2 = getNeighbour(position, mu, REAR);

                        LOOP_DIRAC(alphap) {
                            LOOP_3(ap) {
                                prod1 = -0.5 *
                                        (identityDirac.m[alpha * 4 + alphap] - gamma[mu].m[alpha * 4 + alphap]) *
                                        (link1.m[ELEM_3X3(a, ap)]) *
                                        (*(ELEM_VEC_POSDC(f, positionp1, alphap, ap)));

                                prod2 = -0.5 *
                                        (identityDirac.m[alpha * 4 + alphap] + gamma[mu].m[alpha * 4 + alphap]) *
                                        (link2.m[ELEM_3X3(a, ap)]) *
                                        (*(ELEM_VEC_POSDC(f, positionp2, alphap, ap)));

                                //	Anti-periodic boundary condition
                                (*(ELEM_VEC_POSDC(g, position, alpha, a))) +=
                                    (mu != T_INDX || position.pos[T_INDX] != (lattice_param.n_T - 1)) ? prod1 : -prod1;
                                (*(ELEM_VEC_POSDC(g, position, alpha, a))) +=
                                    (mu != T_INDX || position.pos[T_INDX] != 0) ? prod2 : -prod2;
                            }
                        }
                    }
                    if (c_SW_dirac_operator != 0.0) {
                        LOOP_DIRAC(alphap) {
                            LOOP_3(ap) {
                                (*(ELEM_VEC_POSDC(g, position, alpha, a))) += (c_SW_dirac_operator / 2.0) * ((*(ELEM_VEC_POS(sigmamunuFmunu_dirac_operator, position))).m[ELEM_DCxDC(alpha, a, alphap, ap)]) * (*(ELEM_VEC_POSDC(f, position, alphap, ap)));
                            }
                        }
                    }
                }
            }
        }
    }
}

double invertDiracOperator(const double kappa, Mtrx3x3 *U, Scalar *source, Scalar *inverse_column, const double tolerance, double (*inversion_algorithm)(void (*)(Scalar *, Scalar *), Scalar *, Scalar *, double, size_t)) {
    //  TODO: CREATE INITIALIZER AND DESTRUCTOR OF DIRAC OPERATOR
    U_dirac_operator = U;
    kappa_dirac_operator = kappa;

    size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    return inversion_algorithm(DiracOperator, source, inverse_column, tolerance, sizeof_vector);
}

void initializePauliTerm(Mtrx3x3 *U, const double c_SW, DiracColorMatrix *sigmamunuFmunu) {
    LorentzIdx mu, nu;
    DiracIdx alpha, beta, delta;
    MtrxIdx3 a, b;
    PosVec position;

    DiracMatrix sigma[DIM][DIM];

    Mtrx3x3 Qmunu, Qnumu;

    LOOP_LORENTZ(mu) {
        LOOP_LORENTZ(nu) {
            LOOP_DIRAC(alpha) {
                LOOP_DIRAC(beta) {
                    sigma[mu][nu].m[alpha * 4 + beta] = 0.0;
                    LOOP_DIRAC(delta) {
                        sigma[mu][nu].m[alpha * 4 + beta] +=
                            (gamma[mu].m[alpha * 4 + delta] *
                                 gamma[nu].m[delta * 4 + beta] -
                             gamma[nu].m[alpha * 4 + delta] *
                                 gamma[mu].m[delta * 4 + beta]) /
                            (2.0 * I);
                    }
                }
            }
        }
    }

    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    LOOP_DIRAC(beta) {
                        LOOP_3(b) {
                            (*(ELEM_VEC_POS(sigmamunuFmunu, position))).m[ELEM_DCxDC(alpha, a, beta, b)] = 0.0;
                        }
                    }
                }
            }

            LOOP_LORENTZ(mu) {
                LOOP_LORENTZ(nu) {
                    CloverTerm(U, position, mu, nu, &Qmunu);
                    CloverTerm(U, position, nu, mu, &Qnumu);
                    LOOP_DIRAC(alpha) {
                        LOOP_3(a) {
                            LOOP_DIRAC(beta) {
                                LOOP_3(b) {
                                    (*(ELEM_VEC_POS(sigmamunuFmunu, position))).m[ELEM_DCxDC(alpha, a, beta, b)] += (c_SW / 2.0) *
                                                                                                                    sigma[mu][nu].m[alpha * 4 + beta] * (-I / (8.0)) *
                                                                                                                    (Qmunu.m[ELEM_3X3(a, b)] -
                                                                                                                     Qnumu.m[ELEM_3X3(a, b)]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

double invertImprovedDiracOperator(const double kappa, Mtrx3x3 *U, DiracColorMatrix *sigmamunuFmunu, Scalar *source, Scalar *inverse_column, const double tolerance, double (*inversion_algorithm)(void (*)(Scalar *, Scalar *), Scalar *, Scalar *, double, size_t)) {
    //  TODO: CREATE INITIALIZER AND DESTRUCTOR OF DIRAC OPERATOR
    U_dirac_operator = U;
    kappa_dirac_operator = kappa;

    size_t sizeof_vector = lattice_param.volume * 4 * Nc;
    sigmamunuFmunu_dirac_operator = sigmamunuFmunu;

    double norm = inversion_algorithm(DiracOperator, source, inverse_column, tolerance, sizeof_vector);
    return norm;
}

void printDiracOperator(FILE *file_dirac_op) {
    //	Imprime a Matriz de Dirac explicitamente (Comparação com Mathematica)

    PosVec position1, position2;

    Scalar *f;
    Scalar *g;

    f = (double complex *)calloc(lattice_param.volume * 4 * Nc, sizeof(double complex));
    g = (double complex *)calloc(lattice_param.volume * 4 * Nc, sizeof(double complex));

    DiracIdx alpha, beta;
    MtrxIdx3 a, b;
    LorentzIdx t1, t2;
    LOOP_TEMPORAL(t1) {
        position1.pos[T_INDX] = t1;
        LOOP_SPATIAL(position1) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    (*(ELEM_VEC_POSDC(f, position1, alpha, a))) = 1.0;
                    DiracOperator(f, g);
                    (*(ELEM_VEC_POSDC(f, position1, alpha, a))) = 0.0;
                    LOOP_TEMPORAL(t2) {
                        position2.pos[T_INDX] = t2;
                        LOOP_SPATIAL(position2) {
                            LOOP_DIRAC(beta) {
                                LOOP_3(b) {
                                    fprintf(file_dirac_op, "%.20lf+I(%.20lf)\t",
                                            creal((*(ELEM_VEC_POSDC(g,
                                                                    position2,
                                                                    beta,
                                                                    b)))),
                                            cimag((*(ELEM_VEC_POSDC(g,
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

void printInverse(const char *restrict inverse_filename, Scalar *restrict inverse[4][Nc]) {
    //	Imprime em um arquivo as colunas da matriz inversa do operador de Dirac

    size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    FILE *inverse_file;
    if (!(inverse_file = fopen(inverse_filename, "wb"))) {
        fprintf(stderr,
                "Error: Problem creating file %s for inverse of Dirac matrix.\n", inverse_filename);
        return;
    }

    DiracIdx alpha;
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

DiracMatrix calculateMatrixSlash(double vector[DIM]) {
    DiracMatrix matrix_slashed;

    short unsigned alpha, beta;
    short unsigned mu;
    LOOP_DIRAC(alpha) {
        LOOP_DIRAC(beta) {
            matrix_slashed.m[alpha * 4 + beta] = 0.0;
            LOOP_LORENTZ(mu) {
                matrix_slashed.m[alpha * 4 + beta] +=
                    vector[mu] * gamma[mu].m[alpha * 4 + beta];
            }
        }
    }

    return matrix_slashed;
}

DiracMatrix colortraceFromDiracColor(DiracColorMatrix m_DC) {
    DiracMatrix m_D;
    unsigned short alpha, beta;
    MtrxIdx3 a;
    LOOP_DIRAC(alpha) {
        LOOP_DIRAC(beta) {
            m_D.m[alpha * 4 + beta] = 0.0;
            LOOP_3(a) {
                m_D.m[alpha * 4 + beta] += m_DC.m[ELEM_DCxDC(alpha, a, beta, a)];
            }
        }
    }
    return m_D;
}

Scalar diractraceFromDirac(DiracMatrix m_D) {
    Scalar trace = 0.0;
    unsigned short alpha;
    LOOP_DIRAC(alpha) {
        trace += m_D.m[alpha * 4 + alpha];
    }
    return trace;
}

void aK_from_ap(double ap[4], double aK[4]) {
    LorentzIdx mu;
    LOOP_LORENTZ(mu) {
        aK[mu] = sin(ap[mu]);
    }
}

void aQ_from_ap(double ap[4], double aQ[4]) {
    LorentzIdx mu;
    LOOP_LORENTZ(mu) {
        aQ[mu] = 2.0 * sin(ap[mu] / 2.0);
    }
}

double dotprod(double ap[4], double aq[4]) {
    double dot_prod = 0.0;
    LorentzIdx mu;
    LOOP_LORENTZ(mu) {
        dot_prod += ap[mu] * aq[mu];
    }

    return dot_prod;
}