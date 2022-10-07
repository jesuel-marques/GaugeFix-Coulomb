#ifndef FIELDS_H
#define FIELDS_H

#include <SU3_parameters.h>
#include <types.h>

//  Does the pointer arithmetic to get the correct index in the gauge transformation

#define GET_GT(position)      (((position.t * N_SPC \
                                        + position.k) * N_SPC \
                                                  + position.j) * N_SPC \
                                                            + position.i) 


//  Does the pointer arithmetic to get the correct index in the configuration



#define GET_LINK_U(position, mu) ((((position.t * N_SPC \
                                            + position.k) * N_SPC \
                                                      + position.j) * N_SPC \
                                                                + position.i) * DIM \
                                                                                + mu)


mtrx_3x3 *get_gaugetransf(mtrx_3x3 * restrict G, const pos_vec position);

in_gt_data_type  *get_gaugetransf_in (in_gt_data_type  * restrict G_in,  const pos_vec position);
out_gt_data_type *get_gaugetransf_out(out_gt_data_type * restrict G_out, const pos_vec position);


mtrx_3x3 * get_link(mtrx_3x3 *U, const pos_vec position, const lorentz_idx mu);

in_cfg_data_type  *get_link_in (in_cfg_data_type  * restrict U_in,  const pos_vec position, const lorentz_idx mu);
out_cfg_data_type *get_link_out(out_cfg_data_type * restrict U_out, const pos_vec position, const lorentz_idx mu);

void get_link_matrix(mtrx_3x3 * restrict U, const pos_vec position, const lorentz_idx mu, direction dir, mtrx_3x3 * restrict u);


// void copy_3x3_config(mtrx_3x3 *U, mtrx_3x3 *U_copy);

double average_det(mtrx_3x3 * restrict U);

short SU3_reunitarize_U_G(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G);
// short SU3_reunitarize_U(mtrx_3x3 *restrict U) ;
// short SU3_reunitarize_G(mtrx_3x3 *restrict G);

#endif