#ifndef MISC_H
#define MISC_H

#define OMP_PARALLEL_BASIC "omp parallel for num_threads(NUM_THREADS) schedule(dynamic)"
#define OMP_PARALLEL_FOR  _Pragma(OMP_PARALLEL_BASIC)

//  used to test if allocation was successful
#define TEST_ALLOCATION(a) test_allocation_function(a, __func__) 
#define GREETER()   greeter_function(__FILE__)


void greeter_function(const char * restrict program_name);

short test_allocation_function(const void *pointer,
                               const char *location);

int write_sweeps_to_gaugefix(char * config_filename,
                             int sweeps);

#endif  //MISC_H