
#include <SU3_ops.h>
#include <config_generation.h>
#include <fields.h>
#include <fields_io.h>
#include <geometry.h>
#include <heatbath.h>
#include <plaquettes.h>
#include <ranlux.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#define MAX_FILE_SIZE 1000

int main(int argc, char** argv) {
    char* config_filename = argv[1];
    char complete_config_filename[MAX_FILE_SIZE];

    initGeometry(atoi(argv[2]), atoi(argv[3]));

    double beta = atof(argv[4]);
    unsigned max_configs = atoi(argv[5]);

    unsigned step_to_save = atoi(argv[6]);
    unsigned thermalization = atoi(argv[7]);

    rlxd_init(1, 970223);
    int* ranlux_state;

    FILE* ranlux_state_file;
    char ranlux_state_filename[MAX_FILE_SIZE];
    sprintf(ranlux_state_filename, "%s_%dx%d_beta_%f_%s_start.rlx_state", config_filename, lattice_param.n_SPC, lattice_param.n_T, beta, argv[8]);

    ranlux_state_file = fopen(ranlux_state_filename, "wb");
    ranlux_state = (int*)malloc(rlxd_size() * sizeof(int));
    // fread(ranlux_state, sizeof(int), rlxd_size(), ranlux_state_file);
    // rlxd_reset(ranlux_state);

    if (!validGeometricParametersQ() || argc != 9) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Input config filename, n_SPC, n_T, beta, max configs, steps to save, thermalization, cold/hot "
                "in this order.\n"
                "Check that:\n"
                "n_SPC > 0 and n_T >0\n");
        fprintf(stderr,
                "\n Provided parameters: \n"
                "n_SPC: %d \nn_T: %d \n",
                lattice_param.n_SPC, lattice_param.n_T);
    }

    Mtrx3x3* U = allocate3x3Field(DIM * lattice_param.volume);
    if (U == NULL) {
        fprintf(stderr, "Could not allocate memory for configs %s.\n",
                config_filename);
        return EXIT_FAILURE;
    }

    if (!strcmp(argv[8], "cold")) {
        setFieldToIdentity(U, lattice_param.volume * DIM);
    } else {
        setFieldSU3Random(U, lattice_param.volume * DIM);
    }

    double av_plaq = averagePlaquette(U, "total");

    unsigned configs_saved = 0;
    for (unsigned long sweeps = 0; configs_saved < max_configs; sweeps++) {
        av_plaq += -updateLattice(U, beta, HeatBathSU3) / (2.0 * Nc * pow(lattice_param.n_SPC, 3.0) * lattice_param.n_T * beta);

        printf("%.10lf\n", av_plaq);

        if (sweeps >= thermalization && (sweeps - thermalization) % step_to_save == 0) {
            sprintf(complete_config_filename, "%s_sweep_%lu_%dx%d_beta_%f_%s_start", config_filename, sweeps, lattice_param.n_SPC, lattice_param.n_T, beta, argv[8]);

            if (writeConfig(U, complete_config_filename)) {
                fprintf(stderr, "Config writing failed for config %s.\n", config_filename);
            } else {
                printf("U written for config %s.\n", config_filename);
                configs_saved++;
            }
        }
    }
    free(U);

    rlxd_get(ranlux_state);
    fwrite(ranlux_state, sizeof(int), rlxd_size(), ranlux_state_file);
    fclose(ranlux_state_file);

    return EXIT_SUCCESS;
}