//	srun -p DevQ -N 1 -A nuim01 -t 1:00:00 --pty bash
//	gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -O4 -march=skylake -fopenmp -w
//	mpiicc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -O3 -ipo -xHASWELL -axSKYLAKE,CASCADELAKE,TIGERLAKE -qopt-zmm-usage=high -qopenmp

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>

// #include <mpi.h>
//Since American English is painful on the eyes
// #define MPI_Finalise() MPI_Finalize()
// #include <omp.h>

#include "SU3_parameters.h"			//	Simulation parameters

#include "source/SU3_ops.h"
#include "source/lattice.h"			//	Initialization functions and calculations of
									//	positions and links on the lattice.

#include "source/gauge_fixing.h"	//	Specific functions involved in the gauge-fixing


#include "source/measurement.h"

// const char configs_dir_name[] = "/home/postgrad/jesuel/configs/";
const char configs_dir_name[] = "./configs/";


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
	// const int nconfig = max_configs;
	// //Calculate the number of configs per rank
	// int config_per_rank = nconfig / size;
	// // The for loop divides the work up manually. Instead of using config++ we iterate by the number of configs per rank
	// for (int config = rank + 1; config <= nconfig; config += size) {
	//#pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic) 
	for (unsigned config = 1; config <= max_configs; config ++) {
		int actual_config_nr = 1000+10*(config-1);
	
		in_cfg_data_type * U_in = (in_cfg_data_type *) malloc(Volume * d * sizeof(in_cfg_data_type));
		test_allocation(U_in, "main");

		SU3_load_config(name_configuration_file(actual_config_nr), U_in);
		byte_swap(U_in, sizeof(float), Volume * d * sizeof(in_cfg_data_type));
		
		matrix_3x3_double * U_double = (matrix_3x3_double *) malloc(Volume * d * sizeof(matrix_3x3_double));
		test_allocation(U_double, "main");
		SU3_convert_config_fd(U_in, U_double);
		free(U_in);

		// pos_vec position = assign_position(0, 0, 0, 0);

		// SU3_print_matrix(get_link(U_double, position, 0), "U");
		// getchar();

    	// position = assign_position(Nxyz - 1, Nxyz - 1 , Nxyz - 1 , Nt - 1);
    
		// SU3_print_matrix(get_link(U_double, position, d - 1), "U");
		// getchar();
		
		//  fix the gauge
		SU3_gauge_fix(U_double, actual_config_nr);

		out_cfg_data_type * U_out = (out_cfg_data_type *) malloc(Volume * d * Nc * Nc * sizeof(out_cfg_data_type));
		test_allocation(U_out, "main");

		SU3_convert_config_df(U_double, U_out);
		free(U_double);

		// write the gauge fixed configuration based on template name
		SU3_print_config(name_configuration_file(actual_config_nr), ".GF", U_out);
		free(U_out);	//	Free memory allocated for the configuration.
		
	}
	// MPI_Finalise();
	return 0;
}
