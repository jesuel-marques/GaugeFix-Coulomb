#include <SU3_ops.h>
#include <bicgstab.h>
#include <dirac.h>
#include <fields.h>
#include <fields_io.h>
#include <geometry.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SIZE 1000

int main(int argc, char* argv[]) {
    const char* config_filename = argv[1];

    initGeometry(atoi(argv[2]), atoi(argv[3]));
    double kappa_critical = 0.1392;

    double kappa = atof(argv[4]);
    double tolerance = atof(argv[5]);
    if ((argc != 6 && argc != 5) || !validGeometricParametersQ()) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Input config filename, gt filename, n_SPC, n_T, "
                "tolerance and omega_OR in this order.\n"
                "Check that:\n"
                "n_SPC > 0 and n_T >0\n"
                "tolerance > 0 \n");
        fprintf(stderr,
                "\n Provided parameters: \n"
                "n_SPC: %d \nn_T: %d \n",
                lattice_param.n_SPC, lattice_param.n_T);
        fprintf(stderr, "Initializing gauge-fixing parameters to default instead.\n");
    }

    const Mtrx3x3* U = allocate3x3Field(DIM * lattice_param.volume);
    if (U == NULL) {
        fprintf(stderr, "Could not allocate memory for config in file %s.\n",
                config_filename);
        return EXIT_FAILURE;
    }
    int error;
    if (error = loadConfig(U, config_filename)) {
        fprintf(stderr, "Loading of file %s failed. Error %d.\n", config_filename, error);
        free(U);
        return EXIT_FAILURE;
    } else {
        printf("Config from file %s loaded.\n", config_filename);
    }
    // setFieldToIdentity(U, DIM * lattice_param.volume);

    double am = 1.0 / (2.0 * kappa) - 1.0 / (2.0 * kappa_critical);

    printf("kappa: %lf \t am: %lf\n", kappa, am);
    size_t sizeof_vector = lattice_param.volume * 4 * Nc;
    Scalar* inverse_column[4][3];

    int alpha;
    MtrxIdx3 a;
#pragma omp parallel for collapse(2)
    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            Scalar* source = (Scalar*)calloc(sizeof_vector, sizeof(Scalar));

            inverse_column[alpha][a] = (Scalar*)calloc(sizeof_vector, sizeof(Scalar));
            initializePointSource(assignPosition(0, 0, 0, 0), alpha, a, source);
            BiCGStab(kappa, U, source, inverse_column[alpha][a],
                     tolerance / sizeof_vector);
            free(source);
        }
    }
    free(U);  //	Free memory allocated for the configuration.
    //  SAVE PROPAGATOR TO FILE
    FILE* inverse_column_file;
    char inverse_column_filename[MAX_SIZE];

    sprintf(inverse_column_filename, "%s-Landau.prop", config_filename);

    if ((inverse_column_file = fopen(inverse_column_filename, "wb")) == NULL) {
        fprintf(stderr, "Error: Problem creating file %s for inverse of Dirac matrix.\n",
                inverse_column_filename);
    }
    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            if (fwrite(inverse_column[alpha][a],
                       sizeof(Scalar),
                       sizeof_vector,
                       inverse_column_file) != sizeof_vector) {
                fprintf(stderr, "Error: Problem writing to file %s for inverse of Dirac matrix.\n", inverse_column_filename);
            }
        }
    }
    fclose(inverse_column_file);

    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            free(inverse_column[alpha][a]);
        }
    }

    return EXIT_SUCCESS;
}