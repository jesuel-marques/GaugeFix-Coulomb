//	gcc -o plaquette_check plaquette_check.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/four_potential.c  source/measurement.c source/fields_io.c -lm -O4 -march=skylake -fopenmp

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


#include <SU3_parameters.h>			//	Simulation parameters


#include <lattice.h>			//	Initialization functions and calculations of
									//	positions and links on the lattice.
#include <fields.h>
#include <gauge_fixing.h>	//	Specific functions involved in the gauge-fixing
#include <SU3_ops.h>
#include <fields_io.h>
#include <misc.h>
#include <measurement.h>

const char extension_config_in[]  = ".cfg";
const char extension_config_out[] = "_clmb.cfg";

const char extension_gt_in[] = "_clmb.gt";
const char extension_gt_out[] = "_clmb.gt";

int config_exception_list[] = {-1};

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

	#else
		//OMP_PARALLEL_FOR

		for (unsigned config = 0; config < MAX_CONFIGS; config ++) {


	#endif 
			int actual_config_nr = FIRST_CONFIG + CONFIG_STEP * config;
			
			if(is_in_exception_list(actual_config_nr)) {
				//	configurations which don't actually exist
				printf("Skiping configuration %d for being in the exception list.", actual_config_nr);
				continue;
			}
					
			mtrx_3x3 * U = (mtrx_3x3 *) calloc(VOLUME * DIM, sizeof(mtrx_3x3));
			TEST_ALLOCATION(U);
			
			SU3_load_config(actual_config_nr, U);
			print_matrix_3x3(get_link(U, assign_position(0, 0, 0, 0), 0), "First link", 10);
		
			print_matrix_3x3(get_link(U, assign_position(N_SPC - 1, N_SPC - 1 , N_SPC - 1 , N_T - 1), DIM - 1), "Last link", 10);

			printf("\n e2 before reunitarization and gauge transformation: %5.4E \n", SU3_calculate_e2(U));
			printf("average determinant before reunitarization: %.15lf\n",	average_det(U));
			//  calculate plaquette average
			printf("Spatial plaquette before reunitarization: %.10lf\n", spatial_plaquette_average(U));
			printf("Temporal plaquette before reunitarization: %.10lf\n", temporal_plaquette_average(U));


			mtrx_3x3 * G = (mtrx_3x3 *) calloc(VOLUME , sizeof(mtrx_3x3));
			TEST_ALLOCATION(G);

			SU3_load_gauge_transf(actual_config_nr, G);	

			printf("F before reunitarization: e2: %5.4E\n",SU3_calculate_F(U));
			printf("theta before reunitarization: %5.4E\n", SU3_calculate_theta(U));

			// SU3_reunitarize_U_G(U, G);

			printf("Spatial plaquette before gauge transformation: %.10lf\n", spatial_plaquette_average(U));
			printf("Temporal plaquette before gauge transformation: %.10lf\n", temporal_plaquette_average(U));
			
			printf("F before gt: e2: %5.4E\n",SU3_calculate_F(U));
			printf("theta before gt: %5.4E\n", SU3_calculate_theta(U));

			SU3_global_update_U(U, G);

			printf("F after gt: %5.4E\n",SU3_calculate_F(U));
			printf("theta after gt: %5.4E\n", SU3_calculate_theta(U));


			printf("average determinant: %.15lf\n",	average_det(U));
			printf("\n e2: %5.4E \n", SU3_calculate_e2(U));

			//  calculate plaquette average
			printf("Spatial plaquette after reunitarization:  %.10lf\n", spatial_plaquette_average(U));
			printf("Temporal plaquette after reunitarization:  %.10lf\n", temporal_plaquette_average(U));

			free(G);
			free(U);
					
		}
	#ifdef MPI_CODE
		MPI_Finalise();
	#endif

	return EXIT_SUCCESS;
}
