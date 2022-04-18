#ifndef FOURVECTORFIELD_H
#define FOURVECTORFIELD_H

void SU3_calculate_A(matrix_3x3_double *U, const pos_vec position, const lorentz_index mu, matrix_3x3_double *A);

void SU3_divergence_A(matrix_3x3_double *U, const pos_vec position, matrix_3x3_double *div_A);

#endif