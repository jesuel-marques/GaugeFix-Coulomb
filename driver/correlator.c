#include <SU3_ops.h>
#include <dirac.h>
#include <geometry.h>
#include <spectrum.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char** argv) {
    const char* inverse_filename = argv[1];

    initGeometry(atoi(argv[2]), atoi(argv[3]));

    const char* correlator_filename = argv[4];

    const char* correlator_type = argv[5];

    if ((argc != 6) || !validGeometricParametersQ()) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Inverse filename, n_SPC, n_T, "
                "correlator filename and correlator type in this order.\n"
                "Check that:\n"
                "n_SPC > 0 and n_T > 0\n");
        fprintf(stderr,
                "\n Provided parameters: \n"
                "n_SPC: %d \nn_T: %d \n",
                lattice_param.n_SPC, lattice_param.n_T);
        return EXIT_FAILURE;
    }

    double* correlator = (double*)malloc(lattice_param.n_T * sizeof(double));

    Scalar* inverse[4][Nc];
    size_t sizeof_vector = lattice_param.volume * 4 * Nc;

    DiracIdx alpha;
    MtrxIdx3 a;
    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            inverse[alpha][a] = (Scalar*)malloc(sizeof_vector * sizeof(Scalar));
        }
    }

    FILE* inverse_file = fopen(inverse_filename, "rb");
    LOOP_DIRAC(alpha) {
        LOOP_3(a) {
            if (fread(inverse[alpha][a], sizeof(Scalar), sizeof_vector, inverse_file) != sizeof_vector) {
                fprintf(stderr, "Error reading from file %s for alpha %hd and a %hd.\n", inverse_filename, alpha, a);
                return EXIT_FAILURE;
            };
        }
    }

    calculateCorrelator(inverse, correlator, correlator_type);
    fclose(inverse_file);

    FILE* correlator_file = fopen(correlator_filename, "wb");
    fwrite(correlator, lattice_param.n_T, sizeof(double), correlator_file);
    fclose(correlator_file);

    return EXIT_SUCCESS;
}