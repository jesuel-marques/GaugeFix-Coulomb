#ifndef FOURVECTORFIELD_H
#define FOURVECTORFIELD_H

void SU3_calculate_A(float complex* U, const pos_vec position, const unsigned short mu, float complex* A);

void SU3_divergence_A(float complex* U, const pos_vec position, float complex* div_A);

#endif