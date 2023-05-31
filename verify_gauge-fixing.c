/*
    Main program to gauge-fix configurations to Coulomb-gauge.

    Copyright (C) 2023  Jesuel Marques

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Contact: jesuel.leal@usp.br

 */

//	srun -p DevQ -N 1 -A nuim01 -t 1:00:00 --pty bash
//	gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/four_potential.c source/fields_io.c  -lm -O4 -march=skylake -fopenmp -w
//	mpiicc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/four_potential.c  -lm -O3 -ipo -xHASWELL -axSKYLAKE,CASCADELAKE,TIGERLAKE -qopt-zmm-usage=high -qopenmp -DMPI_CODE

#include <SU3_ops.h>
#include <fields.h>
#include <fields_io.h>
#include <gauge_fixing.h>
#include <geometry.h>
#include <plaquette.h>
#include <mpi.h>
#include <stdio.h>  //	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#define MAX_LENGTH_NAME 1000

int main(int argc, char *argv[]) {
    // Starts MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    // MPI needs a communicator to know how to send/receive data. We aren't sending or receiving things here
    MPI_Comm comm = MPI_COMM_WORLD;
    // The rank is the process number
    MPI_Comm_rank(comm, &rank);
    // The size is the number of processes
    MPI_Comm_size(comm, &size);
    ORGaugeFixingParameters gfix_param;

    initGeometry(atoi(argv[3]), atoi(argv[4]));

    gfix_param = initORGaugeFixingParameters(argv[5]);

    if (!rank) {
        if (!validORGaugeFixingParametersQ(gfix_param) ||
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
    }
    MPI_Barrier(comm);

    const int start_config = atoi(argv[6]);
    const int skip_config = atoi(argv[7]);
    const int nconfig = atoi(argv[8]);
    // Calculate the number of configs per rank
    int config_per_rank = nconfig / size;

    for (int config = rank; config < nconfig; config += size) {
        int actual_config_nr = start_config + config * skip_config;

        char config_filename[MAX_LENGTH_NAME];
        char gauge_transf_filename[MAX_LENGTH_NAME];

        sprintf(config_filename, "%s/Run1_%d.cfg", argv[1], actual_config_nr);
        sprintf(gauge_transf_filename, "%s/Run1_%d_e2_%3.2E.gt", argv[2], actual_config_nr, gfix_param.generic_gf.tolerance);

        Mtrx3x3 *U = allocate3x3Field(DIM * lattice_param.volume);
        if (U == NULL) {
            fprintf(stderr, "Could not allocate memory for config in file %s.\n",
                    config_filename);
        }

        if (loadConfig(U, config_filename)) {
            fprintf(stderr, "Loading of file %s failed.\n", config_filename);
            free(U);
        } else {
            printf("Config from file %s loaded.\n", config_filename);
        }

        //	Reunitarizing right away because of loss of precision due to
        //	storing config in single precision.

        if (reunitarizeField(U, DIM * lattice_param.volume)) {
            fprintf(stderr, "Configuration in file %s could not be reunitarized.\n",
                    config_filename);
            free(U);
        }

        Mtrx3x3 *G = allocate3x3Field(lattice_param.volume);
        if (G == NULL) {
            fprintf(stderr,
                    "Could not allocate memory for gauge transformation "
                    "for file %s.\n",
                    config_filename);
            free(U);
        }

        if (loadGaugeTransf(G, gauge_transf_filename)) {
            fprintf(stderr, "Loading of file %s failed.\n", config_filename);
            free(U);
            free(G);
        } else {
            printf("Gauge-transformation from file %s loaded OK.\n", config_filename);
        }

        if (reunitarizeField(G, lattice_param.volume)) {
            fprintf(stderr, "Gauge-transformation in file %s could not be reunitarized.\n",
                    gauge_transf_filename);
            free(U);
            free(G);
        }

        double plaquettes_before = averagePlaquette(U, "total");
        applyGaugeTransformationU(U, G);
        double plaquettes_after = averagePlaquette(U, "total");

        printf("Difference in plaquettes for config %s: %3.2E \n", config_filename, plaquettes_after - plaquettes_before);

        double residue = gfix_param.generic_gf.gfix_proxy(U);
        printf("%d: %3.2E\n", actual_config_nr, residue);
        if (residue > gfix_param.generic_gf.tolerance)
            fprintf(stderr, "WARNING: CONFIG %s not fixed. residue: %3.2E\n", config_filename, residue);

        free(U);  //	Free memory allocated for the configuration.
        free(G);  //	Free memory allocated for gauge transformation.
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}