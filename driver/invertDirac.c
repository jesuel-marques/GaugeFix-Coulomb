#include <SU3_ops.h>
#include <bicgstab.h>
#include <dirac.h>
#include <fields.h>
#include <fields_io.h>
#include <gauge_fixing.h>
#include <geometry.h>
#include <plaquette.h>
#include <ranlux.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SIZE 1000

// #define KAPPA_CRITICAL 0.1392
#define KAPPA_CRITICAL 0.125

int main(int argc, char* argv[]) {
    const char* config_filename = argv[1];

    initGeometry(atoi(argv[2]), atoi(argv[3]));

    const double kappa = atof(argv[4]);
    const double c_SW = atof(argv[5]);

    const double tolerance = atof(argv[6]);
    const char* inverse_filename = argv[7];

    if ((argc != 8) || !validGeometricParametersQ()) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Input config filename, n_SPC, n_T, "
                "kappa, c_SW, tolerance and inverse filename in this order.\n"
                "Check that:\n"
                "n_SPC > 0 and n_T > 0\n"
                "kappa > 0 and tolerance > 0 \n");
        fprintf(stderr,
                "\n Provided parameters: \n"
                "n_SPC: %d \nn_T: %d \n",
                lattice_param.n_SPC, lattice_param.n_T);
        return EXIT_FAILURE;
    }

    Mtrx3x3* U = allocate3x3Field(DIM * lattice_param.volume);
    if (U == NULL) {
        fprintf(stderr, "Could not allocate memory for config in file %s.\n",
                config_filename);
        return EXIT_FAILURE;
    }

    setFieldToIdentity(U, DIM * lattice_param.volume);

    // int error;
    // if ((error = loadConfig(U, config_filename))) {
    //     fprintf(stderr, "Loading of file %s failed. Error %d.\n", config_filename, error);
    //     free(U);
    //     return EXIT_FAILURE;
    // } else {
    //     printf("Config from file %s loaded.\n", config_filename);
    // }

    // Mtrx3x3* G = allocate3x3Field(lattice_param.volume);
    // if (G == NULL) {
    //     fprintf(stderr, "Could not allocate memory for random gauge transformation.\n");
    //     return EXIT_FAILURE;
    // }
    // printf("Plaquette average before GT: %.18lf. e2: %3.2E \n", creal(averagePlaquette(U, "total")), calculate_e2(U, LANDAU));
    // rlxd_init(1, 12);
    // setFieldSU3Random(G, lattice_param.volume);
    // applyGaugeTransformationU(U, G);
    // free(G);

    // loadConfigPlainText(U, config_filename);
    printf("Plaquette average after GT: %.18lf. e2: %3.2E \n", creal(averagePlaquette(U, "total")), calculate_e2(U, LANDAU));
    // Mtrx3x3 link;
    // PosVec position_to_print = assignPosition(atoi(argv[8]), atoi(argv[9]), atoi(argv[10]), atoi(argv[11]));
    // getLinkMatrix(U, position_to_print, atoi(argv[12]), REAR, &link);
    // printMatrix3x3(&link);
    // getchar();

    double am = 1.0 / (2.0 * kappa) - 1.0 / (2.0 * KAPPA_CRITICAL);
    printf("kappa: %lf \t am: %lf\n", kappa, am);

    size_t sizeof_vector = lattice_param.volume * 4 * Nc;
    Scalar* restrict inverse_columns[4][Nc];

    DiracColorMatrix* pauli_term = NULL;
    if (c_SW != 0.0) {
        pauli_term =
            (DiracColorMatrix*)calloc(sizeof_vector, sizeof(DiracColorMatrix));

        initializePauliTerm(U, c_SW, pauli_term);
    }

    DiracIdx alpha;
    MtrxIdx3 a;
#pragma omp parallel for collapse(2)
    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            Scalar* source = (Scalar*)calloc(sizeof_vector, sizeof(Scalar));

            inverse_columns[alpha][a] = (Scalar*)calloc(sizeof_vector, sizeof(Scalar));
            initializePointSource(assignPosition(0, 0, 0, 0), alpha, a, source);
            if (c_SW == 0.0) {
                invertDiracOperator(kappa, U, source,
                                    inverse_columns[alpha][a],
                                    tolerance * sizeof_vector, BiCGStab);
            } else {
                invertImprovedDiracOperator(kappa, U, pauli_term, source,
                                            inverse_columns[alpha][a],
                                            tolerance * sizeof_vector, BiCGStab);
            }
            free(source);
        }
    }

    FILE* dirac_op_file;
    char dirac_op_filename[MAX_SIZE];
    sprintf(dirac_op_filename, "%s.dirac_op", config_filename);
    dirac_op_file = fopen(dirac_op_filename, "w");
    printDiracOperator(dirac_op_file);
    fclose(dirac_op_file);

    free(U);  //	Free memory allocated for the configuration.
    if (c_SW != 0) {
        free(pauli_term);
    }

    //  Save inverse to file

    printInverse(inverse_filename, inverse_columns);

    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            free(inverse_columns[alpha][a]);
        }
    }

    return EXIT_SUCCESS;
}