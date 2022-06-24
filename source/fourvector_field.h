#ifndef FOURVECTORFIELD_H
#define FOURVECTORFIELD_H

void SU3_calculate_A(mtrx_3x3 * restrict U, const pos_vec position, const lorentz_idx mu, 
                                                                    mtrx_3x3 * restrict A);

void SU3_divergence_A(mtrx_3x3 * restrict U, const pos_vec position, 
                                                                mtrx_3x3 * restrict div_A);

#endif