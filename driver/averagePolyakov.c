
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
#include <geometry.h>
#include <plaquette.h>
#include <polyakov.h>
#include <ranlux.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#define MAX_LENGTH_NAME 100

int main(int argc, char* argv[]) {
    char* config_filename = argv[1];

    initGeometry(atoi(argv[2]), atoi(argv[3]));

    if ((argc != 4) ||
        !validGeometricParametersQ()) {
        fprintf(stderr,
                "Bad input.\n"
                "Usage: Input config filename, n_SPC, n_T \n"
                "Check that:\n"
                "n_SPC > 0 and n_T >0\n");
        fprintf(stderr,
                "\n Provided parameters: \n"
                "n_SPC: %d \nn_T: %d \n",
                lattice_param.n_SPC, lattice_param.n_T);
        fprintf(stderr, "Initializing gauge-fixing parameters to default instead.\n");
        return EXIT_FAILURE;
    }

    Mtrx3x3* U = allocate3x3Field(DIM * lattice_param.volume);
    if (U == NULL) {
        fprintf(stderr, "Could not allocate memory for config in file %s.\n",
                config_filename);
        return EXIT_FAILURE;
    }
    // rlxd_init(1, 11111);

    // setFieldSU3Random(U, DIM * lattice_param.volume);

    if (loadConfig(U, config_filename)) {
        fprintf(stderr, "Loading of file %s failed.\n", config_filename);
        free(U);
        return EXIT_FAILURE;
    } else {
        // printf("Config from file %s loaded.\n", config_filename);
    }

    //	Reunitarizing right away because of loss of precision .

    if (reunitarizeField(U, DIM * lattice_param.volume)) {
        fprintf(stderr, "Configuration in file %s could not be reunitarized.\n",
                config_filename);
        free(U);
        return EXIT_FAILURE;
    }

    // doCenterTransformation(U, 0, TWO_PI_OVER_THREE);

    Scalar average_polyakov_loop = averagePolyakovLoop(U);
    printf("%lf+I*(%lf)\n", creal(average_polyakov_loop), cimag(average_polyakov_loop));
    return EXIT_SUCCESS;
}