#ifndef FIELDS_H
#define FIELDS_H

#include <flags.h>
#include <SU3_parameters.h>
#include <types.h>


//  Does the pointer arithmetic to get the correct index in the gauge transformation
#define GET_GT(position)      (((position.t * n_SPC \
                                        + position.k) * n_SPC \
                                                  + position.j) * n_SPC \
                                                            + position.i) 

//  Does the pointer arithmetic to get the correct index in the configuration
#define GET_LINK_U(position, mu) ((((position.t * n_SPC \
                                            + position.k) * n_SPC \
                                                      + position.j) * n_SPC \
                                                                + position.i) * DIM \
                                                                                + mu)

Mtrx3x3 *allocate_field(int elements, 
                        int size_of_elements);

void set_field_to_identity(Mtrx3x3 * restrict field,
                           int elements);

void copy_field(Mtrx3x3 * restrict field, 
                int elements, 
                Mtrx3x3 * restrict field_copy);

int reunitarize_field(Mtrx3x3 * restrict field,
                      int elements);

double average_field_det(Mtrx3x3 *restrict field,
                         int elements) ;



Mtrx3x3 * get_link(Mtrx3x3 *U, 
                   const PosVec position,
                   const LorentzIdx mu);

void get_link_matrix(Mtrx3x3 * restrict U,
                     const PosVec position,
                     const LorentzIdx mu,
                     Direction dir,
                     Mtrx3x3 * restrict u);

Mtrx3x3 *get_gaugetransf(Mtrx3x3 * restrict G,
                         const PosVec position);

#endif  //FIELDS_H