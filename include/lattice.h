#ifndef LATTICE_H
#define LATTICE_H

#include <SU3_parameters.h>
#include <types.h>


#define POSITION_IS_ODD(position)    (((position.i) ^ (position.j) ^ (position.k) ^ (position.t)) & 1)
#define POSITION_IS_EVEN(position)   !(POSITION_IS_ODD(position))


//  if position is odd, then the XOR of the first bit of each element
//  of position must be 1. Take AND with 1 select this first bit. Take the NOT of 
//  the odd code, because want 1 for even and 0 for odd.

//  Odd position means that the sum of the coordinates is odd and equivalente for even


//  Associations between the numeric indices and the lorentz directions

#define X_INDX 0
#define Y_INDX 1
#define Z_INDX 2
#define T_INDX 3

pos_vec assign_position(const pos_index x, const pos_index y, const pos_index z, const pos_index t);

void print_pos_vec(const pos_vec u);

pos_vec hop_pos_plus (const pos_vec u, const lorentz_idx mu),
        hop_pos_minus(const pos_vec u, const lorentz_idx mu);

#endif