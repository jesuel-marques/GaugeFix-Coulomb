#ifndef CONFIGIO_H
#define CONFIGIO_H

#include <stdlib.h>
#include "lattice.h"

#include "../SU3_parameters.h"

void greeter_function(char * program_name);

void handle_input(int argc, char *argv[]);

void create_output_directory(void);

bool is_in_exception_list(const int config_nr);

void SU3_load_config(const unsigned config_nr, mtrx_3x3 *U),
     SU3_write_config(const unsigned config_nr, mtrx_3x3 *U);

void SU3_load_gauge_transf(const unsigned config_nr, mtrx_3x3 * restrict G),
          SU3_write_gauge_transf(const unsigned config_nr, mtrx_3x3 *G);


#endif