
#include <SU3_ops.h>
#include <config_generation.h>
#include <fields.h>
#include <fields_io.h>
#include <geometry.h>
#include <heatbath.h>
#include <plaquette.h>
#include <polyakov.h>
#include <ranlux.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#define MAX_FILE_SIZE 1000
#define MAX_SWEEPS 10000000

int main(int argc, char** argv) {
    char* config_filename = argv[1];
    char complete_config_filename[MAX_FILE_SIZE];

    initGeometry(atoi(argv[2]), atoi(argv[3]));

    double beta = atof(argv[4]);
    unsigned max_configs = atoi(argv[5]);

    unsigned step_to_save = atoi(argv[6]);
    unsigned thermalization = atoi(argv[7]);

    FILE* ranlux_state_file;
    char ranlux_state_filename[MAX_FILE_SIZE];
    sprintf(ranlux_state_filename, "%s_%dx%d_beta_%f_%s_start.rlx_state", config_filename, lattice_param.n_SPC, lattice_param.n_T, beta, argv[8]);

    rlxd_init(1, 32244000);
    int* ranlux_state;
    ranlux_state = (int*)malloc(rlxd_size() * sizeof(int));

    if (!validGeometricParametersQ() || argc < 9) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Input config filename, n_SPC, n_T, beta, max configs, steps to save, thermalization, cold/hot/restart, (initial sweep number) "
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

    unsigned long configs_to_generate = (argc == 10 ? atoi(argv[9]) : 0);

    if (!strcmp(argv[8], "cold")) {
        setFieldToIdentity(U, lattice_param.volume * DIM);
    } else if (!strcmp(argv[8], "hot")) {
        setFieldSU3Random(U, lattice_param.volume * DIM);
    } else if (!strcmp(argv[8], "restart")) {
        sprintf(complete_config_filename, "%s_sweep_%lu_%dx%d_beta_%f_cold_start", config_filename, sweeps, lattice_param.n_SPC, lattice_param.n_T, beta);

        printf("Restarting from config %s.\n", complete_config_filename);
        if (loadConfig(U, complete_config_filename)) {
            fprintf(stderr, "Config reading failed for config %s.\n", complete_config_filename);
            return EXIT_FAILURE;
        }
        ranlux_state_file = fopen(ranlux_state_filename, "rb");

        if (fread(ranlux_state, sizeof(int), rlxd_size(), ranlux_state_file) != rlxd_size()) {
            fprintf(stderr, "Could not read ranlux state from file %s.\n", ranlux_state_filename);
            fclose(ranlux_state_file);
            return EXIT_FAILURE;
        }
        fclose(ranlux_state_file);
        rlxd_reset(ranlux_state);
        sweeps++;
    } else {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Input config filename, n_SPC, n_T, beta, max configs, steps to save, thermalization, cold/hot/restart, (initial sweep number) "
                "in this order.\n"
                "Check that:\n"
                "n_SPC > 0 and n_T >0\n");
        fprintf(stderr,
                "\n Provided parameters: \n"
                "n_SPC: %d \nn_T: %d \n",
                lattice_param.n_SPC, lattice_param.n_T);
        return EXIT_FAILURE;
    }

    double av_plaq = averagePlaquette(U, "total");
    Scalar average_polyakov_loop = averagePolyakovLoop(U);

    unsigned configs_saved = 0;

    printf("sweep: %ld \t plaq: %.10lf \t polyakov loop: %.10lf+I*(%.10lf)\n", sweeps, av_plaq, creal(average_polyakov_loop), cimag(average_polyakov_loop));

    for (; configs_saved < max_configs && sweeps < MAX_SWEEPS; sweeps++) {
        av_plaq = averagePlaquette(U, "total");
        average_polyakov_loop = averagePolyakovLoop(U);

                printf("sweep: %ld \t plaq: %.10lf \t polyakov loop: %.10lf+I*(%.10lf)\n", sweeps, av_plaq, creal(average_polyakov_loop), cimag(average_polyakov_loop));

        sprintf(complete_config_filename, "%s_sweep_%lu_%dx%d_beta_%f_%s_start", config_filename, sweeps, lattice_param.n_SPC, lattice_param.n_T, beta, argv[8]);

        if (writeConfig(U, complete_config_filename)) {
            fprintf(stderr, "Config writing failed for config %s.\n", complete_config_filename);
        } else {
            printf("U written for config %s.\n", config_filename);
            configs_saved++;

            ranlux_state_file = fopen(ranlux_state_filename, "wb");
            rlxd_get(ranlux_state);
            fwrite(ranlux_state, sizeof(int), rlxd_size(), ranlux_state_file);
            fclose(ranlux_state_file);
        }
    }
    free(U);
    free(ranlux_state);

    return EXIT_SUCCESS;
}