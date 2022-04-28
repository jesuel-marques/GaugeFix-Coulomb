//	gcc -o plaquette_check plaquette_check.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  source/measurement.c -lm -O4 -march=skylake -fopenmp

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>

#ifdef MPI_CODE

#include <mpi.h>
//Since American English is painful on the eyes
#define MPI_Finalise() MPI_Finalize()

#endif


#include "SU3_parameters.h"			//	Simulation parameters


#include "source/lattice.h"			//	Initialization functions and calculations of
									//	positions and links on the lattice.

#include "source/gauge_fixing.h"	//	Specific functions involved in the gauge-fixing
#include "source/SU3_ops.h"

#include "source/measurement.h"


const char extension_in[] = "_clmbgf.cfg";
// char extension_in[]  = ".cfg";
char extension_out[] = "_clmbgf.cfg";

const char configs_dir_name[MAX_LENGTH_NAME];	//	input from command line
const char config_template[MAX_LENGTH_NAME] ;	//	input from command line

int main(int argc, char *argv[]) {

	handle_input(argc, argv);

	#ifdef MPI_CODE
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
	// The for loop divides the work up manually. Instead of using config++ we iterate by the number of configs per rank
	for (int config = rank + 1; config <= nconfig; config += size) {
	#pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic)

	#else

		for (unsigned config = 0; config < MAX_CONFIGS; config ++) {

	#endif 
			int actual_config_nr = FIRST_CONFIG + CONFIG_STEP * config;
			
			if(is_in_exception_list(actual_config_nr)) {
				//	configurations which don't actually exist
				printf("Skiping configuration %d for being in the exception list.", actual_config_nr);
				continue;
			}
			
			mtrx_3x3 * U = (mtrx_3x3 *) malloc(VOLUME * DIM * sizeof(mtrx_3x3));
			TEST_ALLOCATION(U);
			
			SU3_load_config(actual_config_nr, U);

			print_matrix_3x3(get_link(U, assign_position(0, 0, 0, 0), 0), "First link", 10);
			getchar();
		
			print_matrix_3x3(get_link(U, assign_position(N_SPC - 1, N_SPC - 1 , N_SPC - 1 , N_T - 1), DIM - 1), "Last link", 10);
			getchar();
			
			//  calculate plaquette average
			printf("Spatial plaquete  %.10lf\n", spatial_plaquette_average(U));
			printf("Temporal plaquete %.10lf\n", temporal_plaquette_average(U));

			check_det_1(U);

			free(U);
		
		}
	#ifdef MPI_CODE
		MPI_Finalise();
	#endif

	return EXIT_SUCCESS;
}
