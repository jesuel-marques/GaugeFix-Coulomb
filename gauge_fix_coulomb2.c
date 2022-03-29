// gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -fopenmp -O4 -march=skylake
// icc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -fast -xCASCADELAKE -g -qopenmp

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include <omp.h>

#include "SU3_parameters.h"			//	Simulation parameters


#include "source/lattice.h"			//	Initialization functions and calculations of
									//	positions and links on the lattice.

#include "source/gauge_fixing.h"	//	Specific functions involved in the gauge-fixing


int main(){

	#pragma omp parallel for num_threads(2) schedule (dynamic) 
		for (int config = 1; config <= 2; config++ ) {
			
			double complex * U = (double complex *) malloc(Volume * d * 3 * 3 * sizeof(double complex));

			if ( U == NULL ) {
				//	Test if allocation was successful.
				printf("Memory allocation failed for configuration or gauge-transformation");
				exit(1); 
			}

			//  load each configuration based on name template
			char config_filename[max_length_name];
			sprintf(config_filename, "configs/NewFormConfig_%d_beta_5.700_Nxyz_%d_Nt_%d.txt",config,Nxyz,Nt);
			//sprintf(config_filename, "configs/Config_%d_beta_6.000_Nxyz_%d_Nt_%d.txt",config,Nxyz,Nt);
			SU3_load_config(config_filename, U);
	
			//  fix the gauge
			SU3_gauge_fix(U, config);

			// print the gauge fixed configuration based on template name
			// print_gaugefixedconfig();

			free(U);	//	Free allocated memory for the configuration.
		
		}
	return 0;
}