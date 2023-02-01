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

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include <mpi.h>

#include <fields.h>
#include <fields_io.h>
#include <gauge_fixing.h>	
#include <geometry.h>
#include <measurement.h>
#include <SU3_ops.h>

#define MAX_LENGTH_NAME 1000

int writeSweepsToGaugefix(char * identifier,
                          int sweeps) {
	
	if(!validGeometricParametersQ()) {
		fprintf(stderr, "Error in geometric parameters\n");
    }

	char filename_sweeps_to_gaugefix[MAX_LENGTH_NAME];
	sprintf(filename_sweeps_to_gaugefix, "sweeps_to_gaugefix_%dx%d.txt", 
									     lattice_param.n_SPC, 
									     lattice_param.n_T);
    FILE* sweeps_to_gaugefix;

    if((sweeps_to_gaugefix = fopen(filename_sweeps_to_gaugefix, "a+")) == NULL) {
        fprintf(stderr, "Error opening file %s to record sweeps needed to gaugefix.\n", 
                        filename_sweeps_to_gaugefix);
        return -1;
    }
    else{
        fprintf(sweeps_to_gaugefix, "%s\t%d\n", identifier, sweeps);
        fclose(sweeps_to_gaugefix);
    }
    return 0;
}

int main(int argc, char *argv[]) {

	//Starts MPI
	MPI_Init(&argc,&argv);
	int rank, size;
	//MPI needs a communicator to know how to send/receive data. We aren't sending or receiving things here
	MPI_Comm comm = MPI_COMM_WORLD;
	//The rank is the process number
	MPI_Comm_rank(comm, &rank);
	//The size is the number of processes
	MPI_Comm_size(comm, &size);
	
	ORGaugeFixingParameters gfix_param;

	initGeometry(atoi(argv[3]), atoi(argv[4]));

	gfix_param = initORGaugeFixingParameters(argv[5]);
		
	

	if(!rank){
		printf("\n%s %s %s %s %s %s\n",argv[1],argv[2],argv[3],argv[4],argv[5], argv[6]);
		if(!validORGaugeFixingParametersQ(gfix_param) ||
			!validGeometricParametersQ()				 ) {
			fprintf(stderr, "Bad input.\n"
							"Usage: Input config filename, gt filename, n_SPC, n_T, "
							"tolerance and omega_OR in this order.\n"
							"Check that:\n"
							"1 < omega_OR < 2\n"
							"n_SPC > 0 and n_T >0\n"
							"tolerance > 0 \n");
			fprintf(stderr, "\n Provided parameters: \n"
							"n_SPC: %d \nn_T: %d \n"
							"omega_OR : %lf \n"
							"tolerance : %3.2E \n", 
							lattice_param.n_SPC, lattice_param.n_T,
							gfix_param.omega_OR, gfix_param.generic_gf.tolerance);
			fprintf(stderr, "Initializing gauge-fixing parameters to default instead.\n");
			gfix_param = initParametersORDefault();
		}

		printf("\nGauge-fixing parameters provided:\n");
		printORGaugeFixingParameters(gfix_param);
		
		
	

	}
		MPI_Barrier(comm);


	const int start_config = atoi(argv[6]);
	const int skip_config = atoi(argv[7]);
	const int nconfig = atoi(argv[8]);
		//Calculate the number of configs per rank
	int config_per_rank = nconfig / size;


	for (int config = rank ; config < nconfig; config += size) {

		int actual_config_nr = start_config + config * skip_config;


		char config_filename[MAX_LENGTH_NAME];
		char gauge_transf_filename[MAX_LENGTH_NAME];
		
		sprintf(config_filename, "%s/Run1_%d.cfg", argv[1], actual_config_nr);
		sprintf(config_filename, "%s/Run1_%d_e2_%3.2E.gt", argv[2], actual_config_nr, gfix_param.generic_gf.tolerance );

		Mtrx3x3 * U = allocate3x3Field(DIM * lattice_param.volume);
		if(U == NULL) {
			fprintf(stderr, "Could not allocate memory for config in file %s.\n",
							config_filename);
		}

		if(loadConfig(U, config_filename)) {
			fprintf(stderr, "Loading of file %s failed.\n", config_filename);
			free(U);
		}
		else{
			printf("Config from file %s loaded.\n", config_filename);
		}

		//	Reunitarizing right away because of loss of precision due to
		//	storing config in single precision.	

		if(reunitarizeField(U, DIM * lattice_param.volume)) {
			fprintf(stderr, "Configuration in file %s could not be reunitarized.\n",
							config_filename);
			free(U);
		}

		Mtrx3x3 * G = allocate3x3Field(lattice_param.volume);
		if(G == NULL) {
			fprintf(stderr, "Could not allocate memory for gauge transformation " 
							"for file %s.\n",
							config_filename);
			free(U);			
		}

		//  fix the gauge
		int sweeps = gaugefixOverrelaxation(U, G, gfix_param);
		
		if(sweeps < 0) {
			free(U);
			free(G);
			switch(sweeps) {
				case -1:
					fprintf(stderr, "Configuration in file %s could not be gauge-fixed "
									"within the maximum number of sweeps passed.\n"	
									"SOR algorithm seems not to work or be slower than "
									"what the user expected for this particular config. \n", 
									config_filename);
					break;

				case -2:
					fprintf(stderr, "Error in parameters for gauge fixing\n");
					break;

				default:
					break;
			}
			
		}

		//	Record the effort to gauge-fix
		if(sweeps >= 0) {
			printf("Sweeps needed to gauge-fix config from file %s: %d. residue: %3.2E \n", 
					config_filename,
					sweeps,
					gfix_param.generic_gf.gfix_proxy(U));
			writeSweepsToGaugefix(config_filename, sweeps);
		}

		free(U);		//	Free memory allocated for the configuration.

		// write the gauge transformation to file
		if(writeGaugeTransf(G, gauge_transf_filename)) {
			fprintf(stderr, "Gauge transformation writing to file %s"
							"failed for configuration %s .\n",
							gauge_transf_filename,
							config_filename);
		}
		else{
			printf("G for config %s written to file %s.\n", 
					config_filename,
					gauge_transf_filename);
		}
		
		free(G);		//	Free memory allocated for gauge transformation.
	}
	MPI_Finalize();
	return EXIT_SUCCESS;
}