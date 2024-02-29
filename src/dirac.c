#include <SU3_ops.h>
#include <bicgstab.h>
#include <dirac.h>
#include <fields.h>
#include <plaquette.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

static Mtrx3x3 *U_dirac_operator = NULL;
static DiracColorMatrix *pauli_term_dirac_operator = NULL;

static double kappa_dirac_operator;
static double c_SW_dirac_operator;

static DiracMatrix *gamma = NULL;

static DiracMatrix *identity_plus_gamma = NULL;
static DiracMatrix *identity_minus_gamma = NULL;

static DiracMatrix *sigma = NULL;

//  static DiracMatrix identity_minus_gamma[DIM] = {{{1, 0, 1, 0,
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

static DiracMatrix identityDirac;

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

void printDiracColorMatrix(DiracColorMatrix u) {
    unsigned short a, b;
    printf("{");
    LOOP_DC(a) {
        printf("{");
        LOOP_DC(b) {
            printf("%.*lf + I*(%.*lf)", DIGITS_MATRIX, creal(u.m[a * 4 * 3 + b]),
                   DIGITS_MATRIX, cimag(u.m[a * 4 * 3 + b]));

            b != 12 - 1 ? printf(", ") : 0;
        }
        a != 12 - 1 ? printf("},\n") : 0;
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

void initializeDiracMatrices(void) {
    DiracMatrix gamma_T = {{0, 0, 1, 0,
                            0, 0, 0, 1,
                            1, 0, 0, 0,
                            0, 1, 0, 0}},
                gamma_Z = {{0, 0, -I, 0,
                            0, 0, 0, I,
                            I, 0, 0, 0,
                            0, -I, 0, 0}},
                gamma_Y = {{0, 0, 0, -1,
                            0, 0, 1, 0,
                            0, 1, 0, 0,
                            -1, 0, 0, 0}},
                gamma_X = {{0, 0, 0, -I,
                            0, 0, -I, 0,
                            0, I, 0, 0,
                            I, 0, 0, 0}};

    gamma = (DiracMatrix *)calloc(DIM, sizeof(DiracMatrix));
    identity_plus_gamma = (DiracMatrix *)calloc(DIM, sizeof(DiracMatrix));
    identity_minus_gamma = (DiracMatrix *)calloc(DIM, sizeof(DiracMatrix));

    sigma = (DiracMatrix *)calloc(DIM * DIM, sizeof(DiracMatrix));

    LorentzIdx mu, nu;
    DiracIdx alpha, beta, delta;

    // Initialize gamma matrices

    gamma[T_INDX] = gamma_T;
    gamma[Z_INDX] = gamma_Z;
    gamma[Y_INDX] = gamma_Y;
    gamma[X_INDX] = gamma_X;

    LOOP_DIRAC(alpha) {
        LOOP_DIRAC(beta) {
            identityDirac.m[alpha * 4 + beta] = (alpha == beta) ? 1.0 : 0.0;
        }
    }

    LOOP_LORENTZ(mu) {
        LOOP_DIRAC(alpha) {
            LOOP_DIRAC(beta) {
                identity_plus_gamma[mu].m[alpha * 4 + beta] = identityDirac.m[alpha * 4 + beta] + gamma[mu].m[alpha * 4 + beta];
                identity_minus_gamma[mu].m[alpha * 4 + beta] = identityDirac.m[alpha * 4 + beta] - gamma[mu].m[alpha * 4 + beta];
            }
        }
    }

    // Initialize sigma matrices
    LOOP_LORENTZ(mu) {
        LOOP_LORENTZ(nu) {
            LOOP_DIRAC(alpha) {
                LOOP_DIRAC(beta) {
                    sigma[mu * DIM + nu].m[alpha * 4 + beta] = 0.0;
                    LOOP_DIRAC(delta) {
                        sigma[mu * DIM + nu].m[alpha * 4 + beta] +=
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
}

void initializePauliTerm(void) {
    LorentzIdx mu, nu;
    DiracIdx alpha, beta;
    MtrxIdx3 a, b;
    PosVec position;

    Mtrx3x3 Qmunu, Qnumu;

    // Calculate Pauli term
    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    LOOP_DIRAC(beta) {
                        LOOP_3(b) {
                            (*(ELEM_VEC_POS(pauli_term_dirac_operator, position))).m[ELEM_DCxDC(alpha, a, beta, b)] = 0.0;
                        }
                    }
                }
            }

            LOOP_LORENTZ(mu) {
                LOOP_LORENTZ(nu) {
                    if (mu < nu) {
                        CloverTerm(U_dirac_operator, position, mu, nu, &Qmunu);
                        CloverTerm(U_dirac_operator, position, nu, mu, &Qnumu);
                        LOOP_DIRAC(alpha) {
                            LOOP_3(a) {
                                LOOP_DIRAC(beta) {
                                    LOOP_3(b) {
                                        (*(ELEM_VEC_POS(pauli_term_dirac_operator, position))).m[ELEM_DCxDC(alpha, a, beta, b)] += (c_SW_dirac_operator / 2.0) *
                                                                                                                                   sigma[mu * DIM + nu].m[alpha * 4 + beta] * (-I / (8.0)) *
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
}

void setupDiracOperator(const double kappa, const double c_SW, Mtrx3x3 *U) {
    U_dirac_operator = U;
    kappa_dirac_operator = kappa;
    c_SW_dirac_operator = c_SW;

    initializeDiracMatrices();
    if (c_SW != 0.0) {
        pauli_term_dirac_operator =
            (DiracColorMatrix *)calloc(lattice_param.volume, sizeof(DiracColorMatrix));
        initializePauliTerm();
    }
}

void destroyDiracMatrices(void) {
    free(gamma);
    free(identity_plus_gamma);
    free(identity_minus_gamma);
    free(sigma);
}

void destroyDiracOperator() {
    destroyDiracMatrices();
    if (pauli_term_dirac_operator != NULL) {
        free(pauli_term_dirac_operator);
    }
    U_dirac_operator = NULL;
    kappa_dirac_operator = 0.0;
    c_SW_dirac_operator = 0.0;
}

void NonNullDiracEntries(DiracMatrix M, DiracIdxPair *non_null_entries) {
    DiracIdx alpha, beta;
    DiracIdxPair *non_null_entries_ptr = non_null_entries;

    LOOP_DIRAC(alpha) {
        LOOP_DIRAC(beta) {
            if (M.m[alpha * 4 + beta] != 0.0) {
                *non_null_entries_ptr = (DiracIdxPair){alpha, beta};
                non_null_entries_ptr++;
            }
        }
    }
}

Scalar *DiracOperator(Scalar *f, Scalar *g) {
    //	Calcula g=D.f, ou seja, a matriz de Dirac agindo no vetor f, e coloca o resultado em g.

    PosVec position;
    PosVec positionp1;
    PosVec positionp2;

    Mtrx3x3 link1, link2;

    MtrxIdx3 a, b;
    DiracIdx alpha, beta;
    LorentzIdx mu;

    double antiperiodic_bc1 = 1.0, antiperiodic_bc2 = 1.0;

    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    // Diagonal term: simply copy content of g to f

                    (*(ELEM_VEC_POSDC(g, position, alpha, a))) = (1.0 / (2.0 * kappa_dirac_operator)) * (*(ELEM_VEC_POSDC(f, position, alpha, a)));

                    // Hopping term with forward (positive mu, 1) and backward
                    // (negative mu, 2) neighbors

                    LOOP_LORENTZ(mu) {
                        //	Anti-periodic boundary condition
                        antiperiodic_bc1 = (mu != T_INDX || position.pos[T_INDX] != (lattice_param.n_T - 1)) ? 1.0 : -1.0;
                        antiperiodic_bc2 = (mu != T_INDX || position.pos[T_INDX] != 0) ? 1.0 : -1.0;

                        getLinkMatrix(U_dirac_operator, position, mu, FRONT, &link1);
                        getLinkMatrix(U_dirac_operator, position, mu, REAR, &link2);

                        positionp1 = getNeighbour(position, mu, FRONT);
                        positionp2 = getNeighbour(position, mu, REAR);

                        LOOP_DIRAC(beta) {
                            LOOP_3(b) {
                                (*(ELEM_VEC_POSDC(g, position, alpha, a))) += -0.5 * antiperiodic_bc1 *
                                                                              (identity_minus_gamma[mu].m[alpha * 4 + beta]) * (link1.m[ELEM_3X3(a, b)]) *
                                                                              (*(ELEM_VEC_POSDC(f, positionp1, beta, b)));

                                (*(ELEM_VEC_POSDC(g, position, alpha, a))) += -0.5 * antiperiodic_bc2 *
                                                                              (identity_plus_gamma[mu].m[alpha * 4 + beta]) * (link2.m[ELEM_3X3(a, b)]) *
                                                                              (*(ELEM_VEC_POSDC(f, positionp2, beta, b)));
                            }
                        }
                    }
                    if (pauli_term_dirac_operator != NULL) {
                        LOOP_DIRAC(beta) {
                            LOOP_3(b) {
                                (*(ELEM_VEC_POSDC(g, position, alpha, a))) +=
                                    (*(ELEM_VEC_POS(pauli_term_dirac_operator, position))).m[ELEM_DCxDC(alpha, a, beta, b)] *
                                    (*(ELEM_VEC_POSDC(f, position, beta, b)));
                            }
                        }
                    }
                }
            }
        }
    }

    return g;
}

inline static Scalar *applyHoppingTerm(Scalar *f, Scalar *g) {
    // Hopping term with forward (positive mu, 1) and backward (negative mu, 2) neighbors
    PosVec position, positionp1, positionp2;
    LorentzIdx mu;
    MtrxIdx3 a, b;
    DiracIdx alpha, beta;

    Mtrx3x3 link1, link2;

    double antiperiodicity1, antiperiodicity2;

    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_LORENTZ(mu) {
                //	Anti-periodic boundary condition
                antiperiodicity1 = (mu != T_INDX || position.pos[T_INDX] != (lattice_param.n_T - 1)) ? 1.0 : -1.0;
                antiperiodicity2 = (mu != T_INDX || position.pos[T_INDX] != 0) ? 1.0 : -1.0;

                getLinkMatrix(U_dirac_operator, position, mu, FRONT, &link1);
                getLinkMatrix(U_dirac_operator, position, mu, REAR, &link2);

                positionp1 = getNeighbour(position, mu, FRONT);
                positionp2 = getNeighbour(position, mu, REAR);
                LOOP_DIRAC(alpha) {
                    LOOP_DIRAC(beta) {
                        LOOP_3(a) {
                            LOOP_3(b) {
                                (*(ELEM_VEC_POSDC(g, position, alpha, a))) +=
                                    -0.5 *
                                    antiperiodicity1 *
                                    (identity_minus_gamma[mu].m[alpha * 4 + beta]) *
                                    (link1.m[ELEM_3X3(a, b)]) *
                                    (*(ELEM_VEC_POSDC(f, positionp1, beta, b)));

                                (*(ELEM_VEC_POSDC(g, position, alpha, a))) +=
                                    -0.5 *
                                    antiperiodicity2 *
                                    (identity_plus_gamma[mu].m[alpha * 4 + beta]) *
                                    (link2.m[ELEM_3X3(a, b)]) *
                                    (*(ELEM_VEC_POSDC(f, positionp2, beta, b)));
                            }
                        }
                    }
                }
            }
        }
    }

    return g;
}

inline static Scalar *applyPauliTerm(Scalar *f, Scalar *g) {
    PosVec position;
    MtrxIdx3 a, b;
    DiracIdx alpha, beta;
    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    LOOP_DIRAC(beta) {
                        LOOP_3(b) {
                            (*(ELEM_VEC_POSDC(g, position, alpha, a))) +=
                                (*(ELEM_VEC_POS(pauli_term_dirac_operator, position))).m[ELEM_DCxDC(alpha, a, beta, b)] *
                                (*(ELEM_VEC_POSDC(f, position, beta, b)));
                        }
                    }
                }
            }
        }
    }
    return g;
}

Scalar *newDiracOperator(Scalar *f, Scalar *g) {
    //	Calcula g=D.f, ou seja, a matriz de Dirac agindo no vetor f, e coloca o resultado em g.

    PosVec position;

    MtrxIdx3 a;
    DiracIdx alpha;

    // printf("%lf\n", kappa_dirac_operator);
    // getchar();

    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    // Diagonal term: simply copy content of g to f with a factor

                    (*(ELEM_VEC_POSDC(g, position, alpha, a))) = (1.0 / (2.0 * kappa_dirac_operator)) * (*(ELEM_VEC_POSDC(f, position, alpha, a)));
                }
            }
        }
    }

    applyHoppingTerm(f, g);
    if (pauli_term_dirac_operator != NULL) {
        applyPauliTerm(f, g);
    }

    return g;
}

double invertDiracOperator(Scalar *source, Scalar *inverse_column, const double tolerance, double (*inversion_algorithm)(Scalar *(*)(Scalar *, Scalar *), Scalar *, Scalar *, double, size_t)) {
    size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    return inversion_algorithm(DiracOperator, source, inverse_column, tolerance, sizeof_vector);
}

Scalar *rotateWaveFunction(Scalar cq, Scalar *psi) {
    Scalar *rotated_psi = (Scalar *)calloc(lattice_param.volume * 4 * Nc, sizeof(Scalar));

    PosVec position;
    PosVec positionp1;
    PosVec positionp2;

    Mtrx3x3 link1, link2;

    MtrxIdx3 a, b;
    DiracIdx alpha, beta;
    LorentzIdx mu;

    double antiperiodic_bc1 = 1.0, antiperiodic_bc2 = 1.0;

    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    // Diagonal term: simply copy content of f to g

                    (*(ELEM_VEC_POSDC(rotated_psi, position, alpha, a))) = (*(ELEM_VEC_POSDC(psi, position, alpha, a)));

                    // -Hopping term with forward (positive mu, 1) and backward
                    // (negative mu, 2) neighbors

                    LOOP_LORENTZ(mu) {
                        //	Anti-periodic boundary condition
                        antiperiodic_bc1 = (mu != T_INDX || position.pos[T_INDX] != (lattice_param.n_T - 1)) ? 1.0 : -1.0;
                        antiperiodic_bc2 = (mu != T_INDX || position.pos[T_INDX] != 0) ? 1.0 : -1.0;

                        getLinkMatrix(U_dirac_operator, position, mu, FRONT, &link1);
                        getLinkMatrix(U_dirac_operator, position, mu, REAR, &link2);

                        positionp1 = getNeighbour(position, mu, FRONT);
                        positionp2 = getNeighbour(position, mu, REAR);

                        LOOP_DIRAC(beta) {
                            LOOP_3(b) {
                                (*(ELEM_VEC_POSDC(rotated_psi, position, alpha, a))) +=
                                    -(cq / 2.0) *
                                    antiperiodic_bc1 *
                                    (gamma[mu].m[alpha * 4 + beta]) *
                                    (link1.m[ELEM_3X3(a, b)]) *
                                    (*(ELEM_VEC_POSDC(psi, positionp1, beta, b)));

                                (*(ELEM_VEC_POSDC(rotated_psi, position, alpha, a))) +=
                                    +(cq / 2.0) *
                                    antiperiodic_bc2 *
                                    (gamma[mu].m[alpha * 4 + beta]) *
                                    (link2.m[ELEM_3X3(a, b)]) *
                                    (*(ELEM_VEC_POSDC(psi, positionp2, beta, b)));
                            }
                        }
                    }
                }
            }
        }
    }

    memcpy(psi, rotated_psi, lattice_param.volume * 4 * Nc * sizeof(Scalar));
    free(rotated_psi);
    return psi;
}

Scalar *scaleWaveFunction(Scalar scale, Scalar *psi) {
    PosVec position;
    MtrxIdx3 a;
    DiracIdx alpha;

    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    (*(ELEM_VEC_POSDC(psi, position, alpha, a))) *= scale;
                }
            }
        }
    }

    return psi;
}

void printDiracOperator(FILE *file_dirac_op) {
    //	Imprime a Matriz de Dirac explicitamente (Comparação com Mathematica)

    PosVec position;

    Scalar *f;
    Scalar *g;

    size_t size_of_vector = lattice_param.volume * 4 * Nc;

    f = (Scalar *)calloc(size_of_vector, sizeof(Scalar));
    g = (Scalar *)calloc(size_of_vector, sizeof(Scalar));

    DiracIdx alpha;
    MtrxIdx3 a;
    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            LOOP_DIRAC(alpha) {
                LOOP_3(a) {
                    (*(ELEM_VEC_POSDC(f, position, alpha, a))) = 1.0;
                    DiracOperator(f, g);
                    (*(ELEM_VEC_POSDC(f, position, alpha, a))) = 0.0;
                    fwrite(g, sizeof(Scalar), size_of_vector, file_dirac_op);
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