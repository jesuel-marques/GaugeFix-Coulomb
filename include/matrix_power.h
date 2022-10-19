#ifndef LAPACK_MATRIX_POWER_H
#define LAPACK_MATRIX_POWER_H

#include <types.h>

int matrix_log_3x3  (mtrx_3x3 * restrict a,  
                            mtrx_3x3 * restrict log_of_a);

int matrix_power_3x3(const mtrx_3x3 *a, const scalar power, 
                           mtrx_3x3 *a_to_power);

#endif