//	srun -p DevQ -N 1 -A nuim01 -t 1:00:00 --pty bash
//	gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/four_potential.c source/fields_io.c  -lm -O4 -march=skylake -fopenmp -w
//	mpiicc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/four_potential.c  -lm -O3 -ipo -xHASWELL -axSKYLAKE,CASCADELAKE,TIGERLAKE -qopt-zmm-usage=high -qopenmp -DMPI_CODE

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>

#include <fields.h>
#include <fields_io.h>
#include <gauge_fixing.h>	
#include <lattice.h>
#include <measurement.h>
#include <misc.h>
#include <SU3_parameters.h>			//	Simulation parameters
#include <SU3_ops.h>


short n_SPC;	//   Spatial lattice size
short n_T;		//   Temporal lattice size

int volume;		//	Number of sites in the lattice
int spatial_volume;	

int amount_of_links;
int amount_of_points;
	
int main(int argc, char *argv[]) {

	char config_filename[MAX_LENGTH_NAME] = "";
	char gauge_transf_filename[MAX_LENGTH_NAME] = "";

	strcpy(config_filename		, argv[1]);
	strcpy(gauge_transf_filename, argv[2]);

	n_SPC = atoi(argv[3]);
	n_T   =	atoi(argv[4]);
	double tolerance =	atof(argv[5]);
	float omega_OR = atof(argv[6]);

	int exit_status;

	volume 	       = n_SPC * n_SPC * n_SPC * n_T;
	spatial_volume = n_SPC * n_SPC * n_SPC;	//	Number of sites in the lattice

	amount_of_links = DIM * volume;
	amount_of_points = volume;
	
	
	// GREETER();		

	Mtrx3x3 * U = allocate_field(amount_of_links, sizeof(Mtrx3x3));
	if(U == NULL){
		fprintf(stderr, "Could not allocate memory for config in file %s.\n", config_filename);
		return -1;
	}

	if(SU3_load_config(U, config_filename)){
		fprintf(stderr, "Loading of file %s failed.\n", config_filename);
		free(U);
		return -1;
	}
	else{
		printf("File %s loaded OK.\n", config_filename);
	}

	//	Reunitarizing right away because of loss precision due to
	//	storing config in single precision.	

	if(reunitarize_field(U,  amount_of_links)){
		fprintf(stderr, "Configuration in file %s or its gauge-transformation could not be reunitarized.\n", config_filename);
		free(U);
		return -1;
	}
	Mtrx3x3 * G = allocate_field(amount_of_points, sizeof(Mtrx3x3));
	if(G == NULL){
		free(U);
		return -1;
	}

	set_field_to_identity(G, amount_of_points);
		
	//  fix the gauge

	int sweeps = SU3_gaugefix_overrelaxation(U, G, config_filename, tolerance, omega_OR);

	if(sweeps == -1){
		fprintf(stderr, "Configuration in file %s could not be gauge-fixed.\n \
						 SOR algorithm seems not to work or be too slow\n", config_filename);
		
		free(U);
		free(G);
		return -1;		
	}

	//	Record the effort to gauge-fix
	
	printf("Sweeps needed to gauge-fix config from file %s: %d. e2: %3.2E \n", config_filename, sweeps, SU3_calculate_e2(U));
	write_sweeps_to_gaugefix(config_filename, sweeps);
	
	// write the gauge fixed configuration to file,
	// if(SU3_write_config(actual_config_nr, U)){
	// 	fprintf(stderr, "Config writing failed for config %d.\n", actual_config_nr);
	// }
	// else{
	// 	printf("U written OK for config %d.\n", actual_config_nr);
	// }

	free(U);		//	Free memory allocated for the configuration.

	// write the gauge transformation to file
	if(SU3_write_gauge_transf(G, gauge_transf_filename)){
		fprintf(stderr, "Gauge transformation writing to file %s failed for configuration %s .\n", gauge_transf_filename, config_filename);
	}
	else{
		printf("G written OK for config %s to file %s.\n", config_filename, gauge_transf_filename);
	}
	
	free(G);		//	Free memory allocated for gauge transformation.
	
	return EXIT_SUCCESS;
}