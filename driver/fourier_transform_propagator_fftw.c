#include <SU3_ops.h>
#include <dirac.h>
#include <fftw3.h>
#include <fourier_transform.h>
#include <geometry.h>
#include <spectrum.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#define BIG_NUMBER 1.0E+3

int main(int argc, char** argv) {
    const char* inverse_filename = argv[1];

    initGeometry(atoi(argv[2]), atoi(argv[3]));
    char* traces_filename = argv[4];

    if ((argc != 5) || !validGeometricParametersQ()) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Inverse filename, n_SPC, n_T, "
                "form factor filename in this order.\n"
                "Check that:\n"
                "n_SPC > 0 and n_T > 0\n");
        fprintf(stderr,
                "\n Provided parameters: \n"
                "n_SPC: %d \nn_T: %d \n",
                lattice_param.n_SPC, lattice_param.n_T);
        return EXIT_FAILURE;
    }

    DiracIdx alpha, beta;
    MtrxIdx3 a, b;

    Scalar* inverse_columns[4][Nc];
    size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            inverse_columns[alpha][a] = (Scalar*)malloc(sizeof_vector * sizeof(Scalar));
        }
    }

    FILE* inverse_file = fopen(inverse_filename, "rb");
    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            if (fread(inverse_columns[alpha][a], sizeof(Scalar), sizeof_vector, inverse_file) != sizeof_vector) {
                fprintf(stderr, "Error reading from file %s for alpha %hd and a %hd.\n", inverse_filename, alpha, a);
                return EXIT_FAILURE;
            }
        }
    }
    fclose(inverse_file);

    // DiracMatrix propagator;

    // LOOP_DIRAC(alpha) {
    //     LOOP_DIRAC(beta) {
    //         propagator.m[beta * 4 + alpha] = 0.0;
    //         LOOP_3(a) {
    //             propagator.m[beta * 4 + alpha] +=
    //                 *ELEM_VEC_POSDC(inverse_columns[alpha][a], assignPosition(1, 0, 0, 0), beta, a);
    //         }
    //     }
    // }

    // puts("propagator:\n");
    // printDiracMatrix(propagator);
    // getchar();

    PosVec position;

    fftw_complex* inverse[4 * Nc][4 * Nc];
    fftw_complex* inverse_momentum_space[4 * Nc][4 * Nc];
    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            LOOP_DIRAC(beta) {
                LOOP_3(b) {
                    inverse[alpha * Nc + a][beta * Nc + b] = fftw_alloc_complex(lattice_param.volume);
                    inverse_momentum_space[alpha * Nc + a][beta * Nc + b] = fftw_alloc_complex(lattice_param.volume);
                }
            }
        }
    }

    short unsigned N_T = lattice_param.n_T;
    short unsigned N_S = lattice_param.n_SPC;

    int dimensions[] = {N_T, N_S, N_S, N_S};

    complex double antiperiodic_temporal_BC;

    fftw_plan forward_plan = fftw_plan_dft(4, dimensions, inverse[0][0], inverse_momentum_space[0][0], FFTW_FORWARD, FFTW_MEASURE);

    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            LOOP_TEMPORAL(position.pos[T_INDX]) {
                antiperiodic_temporal_BC = cexp((-2.0 * M_PI * I * (0.5) * position.pos[T_INDX]) / lattice_param.n_T);
                LOOP_SPATIAL(position) {
                    LOOP_DIRAC(beta) {
                        LOOP_3(b) {
                            *(ELEM_VEC_POS(inverse[beta * Nc + b][alpha * Nc + a], position)) = *(ELEM_VEC_POSDC(inverse_columns[alpha][a], position, beta, b)) * antiperiodic_temporal_BC;
                            *(ELEM_VEC_POS(inverse_momentum_space[beta * Nc + b][alpha * Nc + a], position)) = 0.0;
                        }
                    }
                }
            }
            free(inverse_columns[alpha][a]);
        }
    }

    // double complex element;

    // LOOP_DIRAC(alpha) {
    //     LOOP_3(a) {
    //         LOOP_DIRAC(beta) {
    //             LOOP_3(b) {
    //                 element = *(ELEM_VEC_POS(inverse[alpha * Nc + a][beta * Nc + b], assignPosition(1, 0, 0, 0)));
    //                 printf("%lf+I*(%lf) ", creal(element), cimag(element));
    //             }
    //         }
    //         printf("\n");
    //     }
    // }

    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            LOOP_DIRAC(beta) {
                LOOP_3(b) {
                    fftw_execute_dft(forward_plan, inverse[alpha * Nc + a][beta * Nc + b], inverse_momentum_space[alpha * Nc + a][beta * Nc + b]);
                }
            }
        }
    }

    // FILE* propagator_momentum_file = fopen(argv[5], "wb");
    // double complex element;
    // LOOP_DIRAC(alpha) {
    //     LOOP_DIRAC(beta) {
    //         element = 0.0;
    //         LOOP_3(a) {
    //             element += *(ELEM_VEC_POS(inverse_momentum_space[alpha * Nc + a][beta * Nc + a], assignPosition(1, 0, 0, 0)));
    //         }
    //         printf("%.10lf+I*(%.10lf) ", creal(element), cimag(element));
    //     }
    //     printf("\n");
    // }

    // LOOP_DIRAC(alpha) {
    //     LOOP_3(a) {
    //         LOOP_DIRAC(beta) {
    //             LOOP_3(b) {
    //                 fwrite(inverse_momentum_space[alpha * Nc + a][beta * Nc + b], sizeof(fftw_complex), lattice_param.volume, propagator_momentum_file);
    //             }
    //         }
    //     }
    // }
    // fclose(propagator_momentum_file);

    fftw_destroy_plan(forward_plan);

    FILE* traces_file;

    traces_file = fopen(traces_filename, "w");
    initializeDiracMatrices();

#pragma omp parallel shared(inverse_momentum_space, traces_file)
    {
        PosVec n_momentum;
        DiracMatrix propagator_p_space;
        double ap[DIM], aK[DIM];
        Scalar trace_scalar, trace_vector_temporal, trace_vector_spatial;
        double re_trace_vector_spatial, im_trace_vector_spatial;
#pragma omp for collapse(4)
        for (int nt = 0; nt < N_T; nt++) {
            for (int nz = 0; nz < N_S; nz++) {
                for (int ny = 0; ny < N_S; ny++) {
                    for (int nx = 0; nx < N_S; nx++) {
                        // if (!((nx == 0 || nx == N_S / 2) && (ny == 0 || ny == N_S / 2) && (nz == 0 || nz == N_S / 2))) {
                        // if (!(nx == N_S / 2 || ny == N_S / 2 || nz == N_S / 2)) {
                        n_momentum.pos[T_INDX] = nt;
                        n_momentum.pos[Z_INDX] = nz;
                        n_momentum.pos[Y_INDX] = ny;
                        n_momentum.pos[X_INDX] = nx;

                        LOOP_DIRAC(alpha) {
                            LOOP_DIRAC(beta) {
                                propagator_p_space.m[alpha * 4 + beta] = 0.0;

                                LOOP_3(a) {
                                    propagator_p_space.m[alpha * 4 + beta] +=
                                        *(ELEM_VEC_POS(inverse_momentum_space[alpha * Nc + a][beta * Nc + a], n_momentum));
                                }
                            }
                        }

                        trace_scalar = diractraceFromDirac(propagator_p_space) / (4.0 * Nc);

                        ap[T_INDX] = (2.0 * M_PI / N_T) * (nt + 0.5);
                        ap[Z_INDX] = 0.0;
                        ap[Y_INDX] = 0.0;
                        ap[X_INDX] = 0.0;
                        aK_from_ap(ap, aK);

                        trace_vector_temporal =
                            diractraceFromDirac(prodDirac(calculateMatrixSlash(aK),
                                                          propagator_p_space)) /
                            (4.0 * Nc * dotprod(aK, aK) * (-I));

                        ap[T_INDX] = 0.0;
                        ap[Z_INDX] = (2.0 * M_PI / N_S) * nz;
                        ap[Y_INDX] = (2.0 * M_PI / N_S) * ny;
                        ap[X_INDX] = (2.0 * M_PI / N_S) * nx;
                        aK_from_ap(ap, aK);

                        trace_vector_spatial =
                            diractraceFromDirac(prodDirac(calculateMatrixSlash(aK),
                                                          propagator_p_space)) /
                            (4.0 * Nc * dotprod(aK, aK) * (-I));

#pragma omp critical
                        {
                            fprintf(traces_file, "{%d,%d,%d,%d}\t%.10lf+I*(%.10lf)\t%.10lf+I*(%.10lf)\t%.10lf+I*(%.10lf)\n",
                                    nt, nz, ny, nx,
                                    creal(trace_vector_temporal),
                                    cimag(trace_vector_temporal),
                                    fabs(re_trace_vector_spatial = creal(trace_vector_spatial)) < BIG_NUMBER ? re_trace_vector_spatial : NAN,
                                    fabs(im_trace_vector_spatial = cimag(trace_vector_spatial)) < BIG_NUMBER ? im_trace_vector_spatial : NAN,
                                    creal(trace_scalar),
                                    cimag(trace_scalar));
                        }
                    }
                    // }
                }
            }
        }
    }
    destroyDiracMatrices();
    fclose(traces_file);

    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            LOOP_DIRAC(beta) {
                LOOP_3(b) {
                    fftw_free(inverse[alpha * Nc + a][beta * Nc + b]);
                    fftw_free(inverse_momentum_space[alpha * Nc + a][beta * Nc + b]);
                }
            }
        }
    }

    return EXIT_SUCCESS;
}