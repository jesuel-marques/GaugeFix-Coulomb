
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

    unsigned max_configs = atoi(argv[4]);

    float beta = 0.0;

    FILE* ranlux_state_file;
    char ranlux_state_filename[MAX_FILE_SIZE];
    sprintf(ranlux_state_filename, "%s_%dx%d_beta_%f_start.rlx_state", config_filename, lattice_param.n_SPC, lattice_param.n_T, beta);

    rlxd_init(1, 32244000);
    int* ranlux_state;
    ranlux_state = (int*)malloc(rlxd_size() * sizeof(int));

    if (!validGeometricParametersQ() || argc < 5) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Input config filename, n_SPC, n_T, max configs"
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

    long config_number = 0;

    setFieldSU3Random(U, lattice_param.volume * DIM);

    double av_plaq = averagePlaquette(U, "total");
    Scalar average_polyakov_loop = averagePolyakovLoop(U);

    unsigned configs_saved = 0;

    printf("config number: %ld \t plaq: %.10lf \t polyakov loop: %.10lf+I*(%.10lf)\n", config_number, av_plaq, creal(average_polyakov_loop), cimag(average_polyakov_loop));

    for (; configs_saved < max_configs && config_number < MAX_SWEEPS; config_number++) {
        av_plaq = averagePlaquette(U, "total");
        average_polyakov_loop = averagePolyakovLoop(U);

        printf("config number: %ld \t plaq: %.10lf \t polyakov loop: %.10lf+I*(%.10lf)\n", config_number, av_plaq, creal(average_polyakov_loop), cimag(average_polyakov_loop));

        setFieldSU3Random(U, lattice_param.volume * DIM);

        sprintf(complete_config_filename, "%s_sweep_%lu_%dx%d_beta_%f_random", config_filename, config_number, lattice_param.n_SPC, lattice_param.n_T, beta);

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