//	srun -p DevQ -N 1 -A nuim01 -t 1:00:00 --pty bash
//	gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c source/config_io.c  -lm -O4 -march=skylake -fopenmp -w
//	mpiicc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -O3 -ipo -xHASWELL -axSKYLAKE,CASCADELAKE,TIGERLAKE -qopt-zmm-usage=high -qopenmp -DMPI_CODE

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <time.h>

#ifdef MPI_CODE

#include <mpi.h>
//Since American English is painful on the eyes
#define MPI_Finalise() MPI_Finalize()

#endif
//#include <omp.h>

#include "SU3_parameters.h"			//	Simulation parameters

#include "source/SU3_ops.h"
#include "source/config_io.h"
#include "source/lattice.h"			//	Initialization functions and calculations of
									//	positions and links on the lattice.

#include "source/gauge_fixing.h"	//	Specific functions involved in the gauge-fixing
#include "source/integpoly_gauge_fixing.h"
#include "source/config_io.h"

#include "source/measurement.h"

const char extension_config_in[]  = ".cfg";
const char extension_config_out[] = "_clmb.cfg";

const char extension_gt_in [] = "_clmb.gt";
const char extension_gt_out[] = "_clmb.gt";

extern char config_template[];

int config_exception_list[] = {-1};	

//	Configurations to be ignored in the gauge-fixing. -1 indicates the end of the list.

int main(int argc, char *argv[]) {


	#ifdef MPI_CODE
	//	If compiling for running with MPI, then these things have to be defined
	//Starts MPI
	MPI_Init(&argc,&argv);
	int rank, size;
	//MPI needs a communicator to know how to send/receive data. We aren't sending or receiving things here
	MPI_Comm comm = MPI_COMM_WORLD;
	//The rank is the process number
	MPI_Comm_rank(comm, &rank);
	//The size is the number of processes
	MPI_Comm_size(comm, &size);
	const int nconfig = MAX_CONFIGS;
	//Calculate the number of configs per rank
	int config_per_rank = nconfig / size;
	if(!rank){
	#endif
	
		GREETER();
		handle_input(argc, argv);

		// if(create_output_directory()){
		// 	printf("Some error ocurred when creating output directory. Exiting.");
		// 	exit(EXIT_FAILURE);
		// }
	
	#ifdef MPI_CODE	
		}
		MPI_Barrier(comm);
		// The for loop divides the work up manually. Instead of using config++ we iterate by the number of configs per rank
		for (int config = rank ; config < nconfig; config += size) {

	#else
	//	If compiling for not running with MPI, then just use a simple for loop
	// omp_parallel_for
		for (unsigned config = 0; config < MAX_CONFIGS; config ++) {
			
	#endif 
			int actual_config_nr = FIRST_CONFIG + CONFIG_STEP * config;
			
			// if(is_in_exception_list(actual_config_nr)) {
			// 	//	list of configurations to be skipped
			// 	printf("Skiping configuration %d for being in the exception list.\n", actual_config_nr);
			// 	continue;
			// }

			mtrx_3x3 * U = (mtrx_3x3 *) calloc(VOLUME * DIM, sizeof(mtrx_3x3));
			if (TEST_ALLOCATION(U)){
				fprintf(stderr,"Could not allocate memory for config %d. Jumping to the next config.\n", actual_config_nr);
				continue;
			}

			if(SU3_load_config(actual_config_nr, U)){
				fprintf(stderr, "Config %d loading failed.\n", actual_config_nr);
				free(U);
				continue;
			}
			else{
				printf("Config %d loaded OK.\n", actual_config_nr);
			}

			//	Reunitarizing straigh away because of loss precision due to
			//	storing config in single precision.	
			
			mtrx_3x3 * G = (mtrx_3x3 *) calloc(VOLUME , sizeof(mtrx_3x3));
			if(TEST_ALLOCATION(G)){
				fprintf(stderr, "Could not allocate memory for gauge-transformation of config %d.\n", actual_config_nr);
				free(U);
				continue;
			}

			init_gauge_transformation(G);

			if(SU3_reunitarize_U_G(U, G)){
				fprintf(stderr, "Configuration %d or its gauge-transformation could not be reunitarized.\n", actual_config_nr);
				free(U);
				free(G);
				continue;
			}
			
			//  fix the gauge

			integpolyakov_gauge_fix(U, G, actual_config_nr);
			
			int sweeps = SU3_gauge_fix(U, G, actual_config_nr);			
			
			integpolyakov_gauge_fix(U, G, actual_config_nr);

			//	checking if a request to stop has been made
			// if(!system("test -f stop_run")){
			// 	printf("Exiting after request to stop.\n");
			// 	free(U);
			// 	free(G);
			// 	remove("stop_run");
			// 	exit(EXIT_SUCCESS);
			// }

			// write the gauge fixed configuration to file,
			// if(SU3_write_config(actual_config_nr, U)){
			// 	fprintf(stderr, "Config writing failed for config %d.\n", actual_config_nr);
			// }
			// else{
			// 	printf("U written OK for config %d.\n", actual_config_nr);
			// }
			// free(U);		//	Free memory allocated for the configuration.

			// write the gauge transformation to file
			if(SU3_write_gauge_transf(actual_config_nr, G)){
				fprintf(stderr, "Gauge transformation writing failed for config %d.\n", actual_config_nr);
			}
			else{
				printf("G written OK for config %d.\n", actual_config_nr);
			}
			
			free(G);		//	Free memory allocated for gauge transformation.
			
			if(sweeps != -1){
				FILE* sweeps_to_gaugefix;
				char filename_sweeps_to_gaugefix[MAX_LENGTH_NAME];
				sprintf(filename_sweeps_to_gaugefix, "sweeps_to_gaugefix_%s_%dx%d.txt", config_template, N_SPC, N_T);

				if((sweeps_to_gaugefix = fopen(filename_sweeps_to_gaugefix, "a+")) == NULL){
			
					fprintf(stderr, "Error opening file %s for config.\n", filename_sweeps_to_gaugefix);
					continue;
				}
				else{
					// fprintf(sweeps_to_gaugefix, "%d\t%d\n", actual_config_nr, sweeps);
					fflush(sweeps_to_gaugefix);	
					fclose(sweeps_to_gaugefix);
				}
			}
			else{
				printf("Configuration %d could not be gauge-fixed.\n", actual_config_nr);
			}
	#ifndef MPI_CODE
		}
	#else
		}
		MPI_Finalise();
	#endif

	return EXIT_SUCCESS;
}