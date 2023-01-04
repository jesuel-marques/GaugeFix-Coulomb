//	gcc -o plaquette_check plaquette_check.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/four_potential.c  source/measurement.c source/fields_io.c -lm -O4 -march=skylake -fopenmp

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>

#include <fields.h>
#include <fields_io.h>
#include <gauge_fixing.h>
#include <lattice.h>
#include <measurement.h>
#include <misc.h>
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

	int exit_status;

	volume 	       = n_SPC * n_SPC * n_SPC * n_T;
	spatial_volume = n_SPC * n_SPC * n_SPC;	//	Number of sites in the lattice

	amount_of_links = DIM * volume;
	amount_of_points = volume;
	
	
	// GREETER();

	Mtrx3x3 * U = (Mtrx3x3 *) calloc(volume * DIM, sizeof(Mtrx3x3));
	if (TEST_ALLOCATION(U)) {
		fprintf(stderr,"Could not allocate memory for config filename %s",
						config_filename);
		return -1;
	}

	if(loadConfig(U, config_filename)) {
		fprintf(stderr, "Loading of file %s failed.\n", config_filename);
		free(U);
		return -1;
	}
	else{
		printf("File %s loaded OK.\n", config_filename);
	}
	
	// printf("First Link\n");
	// printMatrix3x3(getLink(U, assignPosition(0, 0, 0, 0), 0), "First link", 10);
	// printf("Last link\n");
	// printMatrix3x3(getLink(U, 
	// 							 assignPosition(n_SPC - 1, n_SPC - 1 , n_SPC - 1 , n_T - 1),
	// 							 DIM - 1));

	// printf("\n e2 before reunitarization and gauge transformation: %5.4E \n", calculate_e2(U));
	// printf("average determinant before reunitarization: %.15lf\n",	average_det(U));
	// //  calculate plaquette average
	// printf("Spatial plaquette before reunitarization: %.10lf\n", averageSpatialPlaquette(U));
	// printf("Temporal plaquette before reunitarization: %.10lf\n", averageTemporalPlaquette(U));


	Mtrx3x3 * G = (Mtrx3x3 *) calloc(volume , sizeof(Mtrx3x3));
	TEST_ALLOCATION(G);

	loadGaugeTransf(G, gauge_transf_filename);	

	// printf("F before reunitarization: e2: %5.4E\n",calculateF(U));
	// printf("theta before reunitarization: %5.4E\n", calculateTheta(U));

	reunitarizeUG(U, G);

	// printf("Spatial plaquette before gauge transformation: %.10lf\n", averageSpatialPlaquette(U));
	// printf("Temporal plaquette before gauge transformation: %.10lf\n", averageTemporalPlaquette(U));
	
	// printf("F before gt: e2: %5.4E\n",calculateF(U));
	// printf("theta before gt: %5.4E\n", calculateTheta(U));

	updateGlobalU(U, G);

	// printf("F after gt: %5.4E\n",calculateF(U));
	// printf("theta after gt: %5.4E\n", calculateTheta(U));


	// printf("average determinant: %.15lf\n",	average_det(U));
	{	
		double e2 = calculate_e2(U);
		if(e2 > tolerance)
			printf("\nWARNING e2: %5.4E TOO LARGE for config in file %s with gauge-transformation in file %s \n", calculate_e2(U), config_filename, gauge_transf_filename);

	}
	//  calculate plaquette average
	// printf("Spatial plaquette after reunitarization:  %.10lf\n", averageSpatialPlaquette(U));
	// printf("Temporal plaquette after reunitarization:  %.10lf\n", averageTemporalPlaquette(U));

	free(G);
	free(U);
			
		

	return EXIT_SUCCESS;
}