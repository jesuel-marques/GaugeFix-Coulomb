#ifndef CONFIGIO_H
#define CONFIGIO_H

#include <stdlib.h>
#include "lattice.h"

#include "../SU3_parameters.h"

enum config_io_error {swap_error = 1, load_error = 2, write_error = 4, command_error = 8, naming_error = 16, allocation_error = 32};

void greeter_function(const char * restrict program_name);

void handle_input(int argc, char *argv[]);

int create_output_directory(void);

bool is_in_exception_list(const int config_nr);

int SU3_load_config (const unsigned config_nr, mtrx_3x3 * restrict U),
    SU3_write_config(const unsigned config_nr, mtrx_3x3 * restrict U);

int SU3_load_gauge_transf (const unsigned config_nr, mtrx_3x3 * restrict G),
    SU3_write_gauge_transf(const unsigned config_nr, mtrx_3x3 * restrict G);


#endif