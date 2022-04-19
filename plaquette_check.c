//	gcc -o plaquette_check plaquette_check.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  source/measurement.c -lm -O4 -march=skylake -fopenmp

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>

// #include <mpi.h>
// //Since American English is painful on the eyes
// #define MPI_Finalise() MPI_Finalize()
// #include <omp.h>

#include "SU3_parameters.h"			//	Simulation parameters


#include "source/lattice.h"			//	Initialization functions and calculations of
									//	positions and links on the lattice.

#include "source/gauge_fixing.h"	//	Specific functions involved in the gauge-fixing
#include "source/SU3_ops.h"

#include "source/measurement.h"

const char configs_dir_name[max_length_name];	//	input from command line
const char config_template[max_length_name] ;	//	input from command line

int main(int argc, char *argv[]) {

	handle_input(argc, argv);

	// //Starts MPI
	// MPI_Init(&argc,&argv);
	// int rank, size;
	// //MPI needs a communicator to know how to send/receive data. We aren't sending or receiving things here
	// MPI_Comm comm = MPI_COMM_WORLD;
	// //The rank is the process number
	// MPI_Comm_rank(comm, &rank);
	// //The size is the number of processes
	// MPI_Comm_size(comm, &size);
	// const int nconfig = max_configs;
	// //Calculate the number of configs per rank
	// int config_per_rank = nconfig / size;
	// // The for loop divides the work up manually. Instead of using config++ we iterate by the number of configs per rank
	// for (int config = rank + 1; config <= nconfig; config += size) {
	//#pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic) 
	for (unsigned config = 1; config <= max_configs; config ++) {
		int actual_config_nr = 1000+10*(config-1);

		matrix_3x3_double * U = (matrix_3x3_double *) malloc(Volume * d * sizeof(matrix_3x3_double));
		test_allocation(U, "main");

		SU3_load_config(actual_config_nr, U);

		print_matrix_3x3(get_link(U, assign_position(0, 0, 0, 0), 0), "U");
		getchar();
    
		print_matrix_3x3(get_link(U, assign_position(Nxyz - 1, Nxyz - 1 , Nxyz - 1 , Nt - 1), d - 1), "U");
		getchar();
		
		//  calculate plaquette average
		printf("Spatial plaquete  %.10lf\n", spatial_plaquette_average(U));
    	printf("Temporal plaquete %.10lf\n", temporal_plaquette_average(U));

		free(U);
	
	}
	// MPI_Finalise();
	return 0;
}
