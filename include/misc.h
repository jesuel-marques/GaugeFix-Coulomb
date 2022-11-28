#ifndef MISC_H
#define MISC_H

#include <bool.h>

#include <SU3_parameters.h>

#define OMP_PARALLEL_BASIC "omp parallel for num_threads(NUM_THREADS) schedule(dynamic)"
#define OMP_PARALLEL_FOR  _Pragma(OMP_PARALLEL_BASIC)

//  used to test if allocation was successful
#define TEST_ALLOCATION(a) test_allocation_function(a, __func__) 
#define GREETER()   greeter_function(__FILE__)

typedef struct{

short n_SPC;	//   Spatial lattice size
short n_T;		//   Temporal lattice size

int volume;		//	Number of sites in the lattice
int spatial_volume;	

int amount_of_links;
int amount_of_points;

char config_filename[MAX_LENGTH_NAME];
char gauge_transf_filename[MAX_LENGTH_NAME];

double tolerance;
float omega_OR;
bool pars_error;
} gauge_fixing_parameters;

void greeter_function(const char * restrict program_name);

short test_allocation_function(const void *pointer,
                               const char *location);

int write_sweeps_to_gaugefix(char * config_filename,
                             int sweeps);

#endif  //MISC_H