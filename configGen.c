
#include <SU3_ops.h>
#include <fields.h>
#include <fields_io.h>
#include <geometry.h>
#include <heatbath.h>
#include <measurement.h>
#include <ranlux.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

int main(int argc, char** argv) {
    const char* config_filename = argv[1];

    initGeometry(atoi(argv[2]), atoi(argv[3]));

    double beta = atof(argv[4]);
    int max_sweeps = atoi(argv[5]);

    rlxd_init(1, 970223);

    if (!validGeometricParametersQ()) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Input config filename, n_SPC, n_T "
                "in this order.\n"
                "Check that:\n"
                "n_SPC > 0 and n_T >0\n");
        fprintf(stderr,
                "\n Provided parameters: \n"
                "n_SPC: %d \nn_T: %d \n",
                lattice_param.n_SPC, lattice_param.n_T);
    }

    const Mtrx3x3* U = allocate3x3Field(DIM * lattice_param.volume);
    if (U == NULL) {
        fprintf(stderr, "Could not allocate memory for config in file %s.\n",
                config_filename);
        return EXIT_FAILURE;
    }

    setFieldToIdentity(U, lattice_param.volume * DIM);
    double av_plaq = averagePlaquette(U, "total");
    for (int sweeps = 0; sweeps < max_sweeps; sweeps++) {
        printf("average_plaq: %.10lf\n", av_plaq);
        av_plaq += -updateLattice(U, beta, HeatBathSU3) / (6 * 8 * 8 * 8 * 8 * beta);
    }

    // if (writeConfig(U, strcat(config_filename, "fixed"))) {
    //     fprintf(stderr, "Fixed config writing failed for config %s.\n", config_filename);
    // } else {
    //     printf("U fixed written for config %s.\n", config_filename);
    // }

    free(U);  //

    return EXIT_SUCCESS;
}