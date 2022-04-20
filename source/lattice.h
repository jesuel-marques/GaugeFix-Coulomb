#ifndef LATTICE_H
#define LATTICE_H

#include "../SU3_parameters.h"
#include "SU3_ops.h"

typedef unsigned short pos_index;

typedef struct {
    pos_index i, j, k;
    pos_index t;
} pos_vec;  //	struct for position vectors

#define POSITION_IS_ODD(position)    ((position.i ^ position.j ^ position.k ^ position.t) & 1)
#define POSITION_IS_EVEN(position)   !((position.i ^ position.j ^ position.k ^ position.t) & 1)

//  if position is odd, then the XOR of the first bit of each element
//  of position must be 1. Take AND with 1 select this first bit. Take the NOT of 
//  the odd code, because want 1 for even and 0 for odd.

typedef unsigned short lorentz_idx;

typedef enum {REAR, FRONT} direction; 

#define x_index 0
#define y_index 1
#define z_index 2
#define t_index 3

typedef mtrx_3x3_float in_cfg_data_type;
typedef mtrx_3x3_float out_cfg_data_type;


pos_vec assign_position(const pos_index x, const pos_index y, const pos_index z, const pos_index t);

void print_pos_vec(const pos_vec u);

pos_vec hop_position_positive(const pos_vec u, const lorentz_idx mu),
        hop_position_negative(const pos_vec u, const lorentz_idx mu);

unsigned short position_is_even(const pos_vec position),
               position_is_odd(const pos_vec position);

void test_allocation(const void * pointer, const char * location );

mtrx_3x3_double * get_link(mtrx_3x3_double *U, const pos_vec position, const lorentz_idx mu);

mtrx_3x3_float *get_link_f(mtrx_3x3_float *U, const pos_vec position, const lorentz_idx mu);

void get_link_matrix(mtrx_3x3_double * U, const pos_vec position, const lorentz_idx mu, direction dir, mtrx_3x3_double * u);

void handle_input(int argc, char *argv[]);

void SU3_load_config(const unsigned config_nr, mtrx_3x3_double *U),
     SU3_write_config(const unsigned config_nr, mtrx_3x3_double *U);

void copy_3x3_config(mtrx_3x3_double *U, mtrx_3x3_double *U_copy);

void SU3_convert_config_fd(mtrx_3x3_float *U_float, mtrx_3x3_double *U_double),
     SU3_convert_config_df(mtrx_3x3_double *U_double, mtrx_3x3_float *U_float);

void SU3_reunitarize(mtrx_3x3_double *U);

/*============================JONIVAR'S CODE===============================*/

void block_swap(int *buffer, size_t length);

void block_swap_double(double *buffer, size_t length);

int byte_swap(void* strip, size_t size, size_t length);

#endif
