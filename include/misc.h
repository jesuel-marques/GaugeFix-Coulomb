#ifndef MISC_H
#define MISC_H

#include <stdbool.h>

#define OMP_PARALLEL_BASIC "omp parallel for num_threads(NUM_THREADS) schedule(dynamic)"
#define OMP_PARALLEL_FOR  _Pragma(OMP_PARALLEL_BASIC)
// #define OMP_PARALLEL_FOR_reduction(op, x) _Pragma(OMP_PARALLEL_BASIC ## " reduction(" ## op:x ## ")")


#define TEST_ALLOCATION(a) test_allocation_function(a, __func__) //  used to test if allocation was successful
#define GREETER()   greeter_function(__FILE__)


void greeter_function(const char * restrict program_name);

void handle_input(int argc, char *argv[]);

// int create_output_directory(void);

bool is_in_exception_list(const int config_nr);

short test_allocation_function(const void *pointer, const char *location);

#endif