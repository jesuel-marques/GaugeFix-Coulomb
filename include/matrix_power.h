#ifndef LAPACK_MATRIX_POWER_H
#define LAPACK_MATRIX_POWER_H

#include <types.h>

int matrix_power_3x3(const Mtrx3x3 *a, 
                     const Scalar power, 
                           Mtrx3x3 *a_to_power);

void matrix_log_3x3  (Mtrx3x3 * restrict a,  
                      Mtrx3x3 * restrict log_of_a);

#endif  //LAPACK_MATRIX_POWER_H