//	srun -p DevQ -N 1 -A nuim01 -t 1:00:00 --pty bash
//	gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -O4 -march=skylake -fopenmp
//	mpiicc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -O3 -ipo -xHASWELL -axSKYLAKE,CASCADELAKE,TIGERLAKE -qopt-zmm-usage=high -qopenmp

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>

// #include <mpi.h>
// //Since American English is painful on the eyes
// #define MPI_Finalise() MPI_Finalize()
#include <omp.h>

#include "SU3_parameters.h"			//	Simulation parameters


#include "source/lattice.h"			//	Initialization functions and calculations of
									//	positions and links on the lattice.

#include "source/gauge_fixing.h"	//	Specific functions involved in the gauge-fixing
#include "source/SU3_ops.h"

const char configs_dir_name[] = "configs/";


int main(void){
// int main(int argc, char *argv[]) {

	// //Starts MPI
	// MPI_Init(&argc,&argv);
	// int rank, size;
	// //MPI needs a communicator to know how to send/receive data. We aren't sending or receiving things here
	// MPI_Comm comm = MPI_COMM_WORLD;
	// //The rank is the process number
	// MPI_Comm_rank(comm, &rank);
	// //The size is the number of processes
	// MPI_Comm_size(comm, &size);
	// const int nconfig = 100;
	// //Calculate the number of configs per rank
	// int config_per_rank = nconfig / size;
	// // The for loop divides the work up manually. Instead of using config++ we iterate by the number of configs per rank
	// for (int config = rank + 1; config <= nconfig; config += size) {
	#pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic) 
	for (unsigned short config = 1; config <= max_configs; config ++) {
		
		float complex * U = (float complex *) malloc(Volume * d * 3 * 3 * sizeof(float complex));
		test_allocation(U, "main");
		SU3_load_config(name_configuration_file(config), U);
		
		//  fix the gauge
		SU3_gauge_fix(U, config);

		// write the gauge fixed configuration based on template name
		SU3_print_config(name_configuration_file(config),"GF", U);

		free(U);	//	Free memory allocated for the configuration.
	}
	// MPI_Finalise();
	return 0;
}
