#ifndef FOURVECTORFIELD_H
#define FOURVECTORFIELD_H

void SU3_calculate_A(mtrx_3x3_double *U, const pos_vec position, const lorentz_idx mu, 
                                                                    mtrx_3x3_double *A);

void SU3_divergence_A(mtrx_3x3_double *U, const pos_vec position, 
                                                                mtrx_3x3_double *div_A);

#endif