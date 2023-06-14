#include <SU3_ops.h>
#include <dirac.h>
#include <fourier_transform.h>
#include <geometry.h>
#include <spectrum.h>
#include <stdio.h>
#include <stdlib.h>
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

    DiracColorMatrix* inverse;
    inverse = (DiracColorMatrix*)malloc(lattice_param.volume * sizeof(DiracColorMatrix));

    Scalar* inverse_columns[4][Nc];
    size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    DiracIdx alpha;
    MtrxIdx3 a;
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

    PosVec position;
    DiracIdx beta;
    MtrxIdx3 b;

    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            LOOP_TEMPORAL(position.pos[T_INDX]) {
                LOOP_SPATIAL(position) {
                    LOOP_DIRAC(beta) {
                        LOOP_3(b) {
                            (ELEM_VEC_POS(inverse, position))->m[ELEM_DCxDC(beta, b, alpha, a)] = *ELEM_VEC_POSDC(inverse_columns[alpha][a], position, beta, b);
                        }
                    }
                }
            }
            free(inverse_columns[alpha][a]);
        }
    }

    FILE* traces_file;
    traces_file = fopen(traces_filename, "w");
    short unsigned N_T = lattice_param.n_T;
    short unsigned N_S = lattice_param.n_SPC;

#pragma omp parallel shared(inverse, traces_file)
    {
        DiracMatrix propagator_p_space;
        double ap[DIM], aK[DIM];
        Scalar trace_scalar, trace_vector;
#pragma omp for collapse(4)
        for (int nt = -N_T / 8; nt <= N_T / 8; nt++) {
            for (int nz = -N_S / 8; nz <= N_S / 8; nz++) {
                for (int ny = -N_S / 8; ny <= N_S / 8; ny++) {
                    for (int nx = -N_S / 8; nx <= N_S / 8; nx++) {
                        ap[T_INDX] = (2.0 * M_PI / N_T) * (nt + 0.5);
                        ap[Z_INDX] = (2.0 * M_PI / N_S) * nz;
                        ap[Y_INDX] = (2.0 * M_PI / N_S) * ny;
                        ap[X_INDX] = (2.0 * M_PI / N_S) * nx;

                        aK_from_ap(ap, aK);
                        // printf("{%lf, %lf, %lf, %lf} \n", ap[T_INDX], ap[Z_INDX], ap[Y_INDX], ap[X_INDX]);

                        propagator_p_space = colortraceFromDiracColor(fourierTransform(inverse, ap));
                        // printDiracMatrix(propagator_p_space);
                        // getchar();
                        trace_scalar = diractraceFromDirac(propagator_p_space);
                        trace_vector =
                            diractraceFromDirac(prodDirac(calculateMatrixSlash(aK),
                                                          propagator_p_space));

                        trace_scalar /= (4 * Nc);
                        trace_vector /= (4 * Nc * dotprod(aK, aK) * (-I));

#pragma omp critical
                        {
                            // fwrite(ap, sizeof(double), 4, traces_file);
                            // fwrite(&trace_scalar, sizeof(Scalar), 1, traces_file);
                            // fwrite(&trace_vector, sizeof(Scalar), 1, traces_file);

                            fprintf(traces_file, "{%d,%d,%d,%d}\t %lf+I*(%lf) \t %lf+I*(%lf)\n",
                                    nt, nz, ny, nx,
                                    creal(trace_vector), cimag(trace_vector),
                                    creal(trace_scalar), cimag(trace_scalar));
                        }
                    }
                }
            }
        }
    }
    fclose(traces_file);
    free(inverse);

    return EXIT_SUCCESS;
}