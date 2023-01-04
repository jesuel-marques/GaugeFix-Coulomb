#ifndef LAPACK_MATRIX_POWER_H
#define LAPACK_MATRIX_POWER_H

#include <SU3_ops.h>
#include <types.h>

int powerMtrx3x3(const Mtrx3x3 *a, 
                     const Scalar power, 
                           Mtrx3x3 *a_to_power);

void logMtrx3x3  (Mtrx3x3 * restrict a,  
                      Mtrx3x3 * restrict log_of_a);

#endif  //LAPACK_MATRIX_POWER_H