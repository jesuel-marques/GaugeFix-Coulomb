//	srun -p DevQ -N 1 -A nuim01 -t 1:00:00 --pty bash
//	gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/four_potential.c source/fields_io.c  -lm -O4 -march=skylake -fopenmp -w
//	mpiicc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/four_potential.c  -lm -O3 -ipo -xHASWELL -axSKYLAKE,CASCADELAKE,TIGERLAKE -qopt-zmm-usage=high -qopenmp -DMPI_CODE

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include <fields.h>
#include <fields_io.h>
#include <gauge_fixing.h>	
#include <lattice.h>
#include <measurement.h>
#include <misc.h>
#include <SU3_ops.h>

#define MAX_LENGTH_NAME 100

int writeSweepsToGaugefix(char * identifier,
                          int sweeps) {
	
	if(!validGeometricParametersQ()){
		fprintf(stderr, "Error in geometric parameters\n");
		exit(EXIT_FAILURE);
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

	const char * config_filename 	   = argv[1];
	const char * gauge_transf_filename = argv[2];

	initGeometry(atoi(argv[3]), atoi(argv[4]));
	
	const ORGaugeFixingParameters gfix_param = initORGaugeFixingParameters(atof(argv[5]), 
															  		   atof(argv[6]));
	
	if(argc != 7 || !validORGaugeFixingParametersQ(gfix_param) 
				 || !validGeometricParametersQ()			) {
		fprintf(stderr, "Bad input.\n"
						"Usage: Input config filename, gt filename, n_SPC, n_T, "
						"tolerance and omega_OR in this order.\n"
						"Check that:\n"
						"1 < omega_OR < 2\n"
						"n_SPC > 0 and n_T >0\n"
						"tolerance > 0 \n");
		fprintf(stderr, "n_SPC: %d \nn_T: %d \n"
		                "omega_OR : %lf \n"
						"tolerance : %3.2E \n", 
						lattice_param.n_SPC, lattice_param.n_T,
						gfix_param.omega_OR, gfix_param.tolerance );

		return EXIT_FAILURE;
	}
	
	const Mtrx3x3 * U = allocateField(lattice_param.amount_of_links, sizeof(Mtrx3x3));
	if(U == NULL) {
		fprintf(stderr, "Could not allocate memory for config in file %s.\n",
						config_filename);
		return EXIT_FAILURE;
	}

	if(loadConfig(U, config_filename)) {
		fprintf(stderr, "Loading of file %s failed.\n", config_filename);
		free(U);
		return EXIT_FAILURE;
	}
	else{
		printf("Config from file %s loaded.\n", config_filename);
	}

	//	Reunitarizing right away because of loss of precision due to
	//	storing config in single precision.	

	if(reunitarizeField(U, lattice_param.amount_of_links)) {
		fprintf(stderr, "Configuration in file %s could not be reunitarized.\n",
						 config_filename);
		free(U);
		return EXIT_FAILURE;
	}

	Mtrx3x3 * G = allocateField(lattice_param.amount_of_points, sizeof(Mtrx3x3));
	if(G == NULL) {
		fprintf(stderr, "Could not allocate memory for gauge transformation " 
						"for file %s.\n",
						config_filename);
		free(U);
		return EXIT_FAILURE;
	}
	
	//  fix the gauge

	int sweeps = gaugefixOverrelaxation(U, G, gfix_param);

	if(sweeps < 0) {
		free(U);
		free(G);
		switch(sweeps) {
			case -1:
				fprintf(stderr, "Configuration in file %s could not be gauge-fixed.\n"	
								"SOR algorithm seems not to work or be too slow\n", 
								config_filename);
				break;

			case -2:
				fprintf(stderr, "Error in parameters for gauge fixing\n");
				break;

			default:
				break;
		}
		
		return EXIT_FAILURE;		
	}

	//	Record the effort to gauge-fix
	if( sweeps >= 0) {
		printf("Sweeps needed to gauge-fix config from file %s: %d. e2: %3.2E \n", 
				config_filename,
				sweeps,
				calculate_e2(U));
		writeSweepsToGaugefix(config_filename, sweeps);
	}

	// write the gauge fixed configuration to file,
	// if(writeConfig(U, strcat(config_filename, "fixed"))) {
	// 	fprintf(stderr, "Fixed config writing failed for config %s.\n", config_filename);
	// }
	// else{
	// 	printf("U fixed written for config %s.\n", config_filename);
	// }

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
	finalizeGeometry();
	return EXIT_SUCCESS;
}