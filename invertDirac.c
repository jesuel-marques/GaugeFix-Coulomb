#include <bicgstab.h>
#include <geometry.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    const char* config_filename = argv[1];

    initGeometry(atoi(argv[2]), atoi(argv[3]));

    ORGaugeFixingParameters gfix_param;

    if ((argc != 6 && argc != 5) ||
        !validORGaugeFixingParametersQ(gfix_param) ||
        !validGeometricParametersQ()) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Input config filename, gt filename, n_SPC, n_T, "
                "tolerance and omega_OR in this order.\n"
                "Check that:\n"
                "1 < omega_OR < 2\n"
                "n_SPC > 0 and n_T >0\n"
                "tolerance > 0 \n");
        fprintf(stderr,
                "\n Provided parameters: \n"
                "n_SPC: %d \nn_T: %d \n"
                "omega_OR : %lf \n"
                "tolerance : %3.2E \n",
                lattice_param.n_SPC, lattice_param.n_T,
                gfix_param.omega_OR, gfix_param.generic_gf.tolerance);
        fprintf(stderr, "Initializing gauge-fixing parameters to default instead.\n");
        gfix_param = initParametersORDefault();
    }

    const Mtrx3x3* U = allocate3x3Field(DIM * lattice_param.volume);
    if (U == NULL) {
        fprintf(stderr, "Could not allocate memory for config in file %s.\n",
                config_filename);
        return EXIT_FAILURE;
    }

    if (loadConfig(U, config_filename)) {
        fprintf(stderr, "Loading of file %s failed.\n", config_filename);
        free(U);
        return EXIT_FAILURE;
    } else {
        printf("Config from file %s loaded.\n", config_filename);
    }

    free(U);  //	Free memory allocated for the configuration.

    return EXIT_SUCCESS;
}