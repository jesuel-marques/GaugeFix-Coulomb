//	srun -p DevQ -N 1 -A nuim01 -t 1:00:00 --pty bash
//	gcc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -O4 -march=skylake -fopenmp -w
//	mpiicc -o gauge_fix_coulomb gauge_fix_coulomb.c source/lattice.c source/SU2_ops.c source/SU3_ops.c source/gauge_fixing.c source/fourvector_field.c  -lm -O3 -ipo -xHASWELL -axSKYLAKE,CASCADELAKE,TIGERLAKE -qopt-zmm-usage=high -qopenmp -DMPI_CODE

#include <stdio.h>					//	Standard header files in C
#include <stdlib.h>
#include <string.h>
#include <complex.h>


#ifdef MPI_CODE

#include <mpi.h>
//Since American English is painful on the eyes
#define MPI_Finalise() MPI_Finalize()

#endif
//#include <omp.h>

#include "SU3_parameters.h"			//	Simulation parameters

#include "source/SU3_ops.h"
#include "source/lattice.h"			//	Initialization functions and calculations of
									//	positions and links on the lattice.

#include "source/gauge_fixing.h"	//	Specific functions involved in the gauge-fixing

#include "source/measurement.h"

// const char extension_in[] = "_clmbgf.cfg";
const char extension_in[]  = ".cfg";
const char extension_out[] = "_clmbgf.cfg";

const char configs_dir_name[MAX_LENGTH_NAME];	//	input from command line
const char config_template[MAX_LENGTH_NAME-30] ;	//	input from command line

int config_exception_list[] = {-1};

int main(int argc, char *argv[]) {

	GREETER();
	handle_input(argc, argv);

	
	char filename_sweeps_to_gaugefix[MAX_LENGTH_NAME];
	sprintf(filename_sweeps_to_gaugefix, "sweeps_to_gaugefix_%s_%dx%d.txt", config_template, N_SPC, N_T);
	

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
	// The for loop divides the work up manually. Instead of using config++ we iterate by the number of configs per rank
	for (int config = rank ; config < nconfig; config += size) {

	#else
	//	If compiling for not running with MPI, then just use a simple for loop
	// #pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic)
		for (unsigned config = 0; config < MAX_CONFIGS; config ++) {
			
	#endif 
			int actual_config_nr = FIRST_CONFIG + CONFIG_STEP * config;
			
			if(is_in_exception_list(actual_config_nr)) {
				//	configurations which don't actually exist
				printf("Skiping configuration %d for being in the exception list.\n", actual_config_nr);
				continue;
			}
			
			mtrx_3x3 * U = (mtrx_3x3 *) calloc(VOLUME * DIM, sizeof(mtrx_3x3));
			TEST_ALLOCATION(U);
			
			SU3_load_config(actual_config_nr, U);

			SU3_reunitarize(U);	
			//	Reunitarizing straigh away because of loss precision due to
			//	storing config in single precision.	

			//  fix the gauge
			
			unsigned sweeps = SU3_gauge_fix(U, actual_config_nr);			
			
			FILE* sweeps_to_gaugefix;
    		if((sweeps_to_gaugefix = fopen(filename_sweeps_to_gaugefix, "a+")) == NULL){
        
				fprintf(stderr, "Error creating file %s for config.\n", filename_sweeps_to_gaugefix);
				exit(EXIT_FAILURE);

			}
			fprintf(sweeps_to_gaugefix, "%d\t%u\n", actual_config_nr, sweeps);
			fflush(sweeps_to_gaugefix);	
			fclose(sweeps_to_gaugefix);

			// write the gauge fixed configuration to file
			SU3_write_config(actual_config_nr, U);
			free(U);		//	Free memory allocated for the configuration.

			
			
		}
	#ifdef MPI_CODE
		MPI_Finalise();
	#endif
	return EXIT_SUCCESS;
}
