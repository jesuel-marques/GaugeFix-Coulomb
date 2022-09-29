#ifndef LAPACK_MATRIX_POWER_H
#define LAPACK_MATRIX_POWER_H

#include <lattice.h>

int matrix_power_3x3(const mtrx_3x3 *a, const work_data_type power, 
                           mtrx_3x3 *a_to_power);

#endif