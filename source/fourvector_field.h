#ifndef FOURVECTORFIELD_H
#define FOURVECTORFIELD_H

void SU3_calculate_A(double complex* U, pos_vec position, int mu, double complex* A);

void SU3_divergence_A(double complex* U, pos_vec position, double complex* div_A);

#endif