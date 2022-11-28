#ifndef FIELDS_H
#define FIELDS_H

#include <flags.h>
#include <SU3_parameters.h>
#include <types.h>


//  Does the pointer arithmetic to get the correct index in the gauge transformation
#define GET_GT(G, position)      G + (((position.t * n_SPC \
                                            + position.k) * n_SPC \
                                                + position.j) * n_SPC \
                                                    + position.i) 

//  Does the pointer arithmetic to get the correct index in the configuration
#define GET_LINK(U, position, mu) U + ((((position.t * n_SPC \
                                            + position.k) * n_SPC \
                                                + position.j) * n_SPC \
                                                    + position.i) * DIM \
                                                        + mu)

Mtrx3x3 *allocate_field(int elements, 
                        int size_of_elements);

int set_field_to_identity(const Mtrx3x3 * restrict field,
                           int elements);

int copy_field(const Mtrx3x3 * restrict field, 
                int elements, 
               const Mtrx3x3 * restrict field_copy);

int reunitarize_field(const Mtrx3x3 * restrict field,
                      int elements);

double average_field_det(const Mtrx3x3 *restrict field,
                         int elements) ;

Mtrx3x3 * get_link(const Mtrx3x3 *U, 
                   const PosVec position,
                   const LorentzIdx mu);

void get_link_matrix(const Mtrx3x3 * restrict U,
                     const PosVec position,
                     const LorentzIdx mu,
                     Direction dir,
                     const Mtrx3x3 * restrict u);

Mtrx3x3 *get_gaugetransf(const Mtrx3x3 * restrict G,
                         const PosVec position);

#endif  //FIELDS_H