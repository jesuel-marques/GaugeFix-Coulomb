#ifndef LATTICE_H
#define LATTICE_H
#include <stdbool.h>

#include "../SU3_parameters.h"


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

//  Odd position means that the sum of the coordinates is odd and equivalente for even

typedef unsigned short lorentz_idx;

typedef enum {REAR, FRONT} direction;

//  Associations between the numeric indices and the lorentz directions

#define X_INDX 0
#define Y_INDX 1
#define Z_INDX 2
#define T_INDX 3

//  Does the pointer arithmetic to get the correct index in the configuration

#define GET_LINK(position, mu) ((((position.t * N_SPC \
                                        + position.k) * N_SPC \
                                            + position.j) * N_SPC \
                                                + position.i) * DIM \
                                                    + mu)

//  Data types definitions

typedef complex float in_data_type; 
typedef complex float out_data_type;
typedef complex double work_data_type;

typedef struct {
   complex double m[Nc * Nc];
} mtrx_3x3_double; 

typedef struct {
    complex float m[Nc * Nc];
} mtrx_3x3_float;

typedef mtrx_3x3_float in_cfg_data_type;
typedef mtrx_3x3_float out_cfg_data_type;
typedef mtrx_3x3_double work_cfg_data_type;

typedef work_cfg_data_type mtrx_3x3;


#define TEST_ALLOCATION(a) test_allocation_in_function(a, __func__) //  used to test if allocation was successful
#define GREETER()   greeter_function(__FILE__)

pos_vec assign_position(const pos_index x, const pos_index y, const pos_index z, const pos_index t);

void print_pos_vec(const pos_vec u);

pos_vec hop_position_positive(const pos_vec u, const lorentz_idx mu),
        hop_position_negative(const pos_vec u, const lorentz_idx mu);


void test_allocation_in_function(const void * pointer, const char * location );


mtrx_3x3 * get_link(mtrx_3x3 *U, const pos_vec position, const lorentz_idx mu);

in_cfg_data_type *get_link_in(in_cfg_data_type *U_in, const pos_vec position, const lorentz_idx mu);
out_cfg_data_type *get_link_out(out_cfg_data_type *U_out, const pos_vec position, const lorentz_idx mu);

void get_link_matrix(mtrx_3x3 * U, const pos_vec position, const lorentz_idx mu, direction dir, mtrx_3x3 * u);

void greeter_function(char * program_name);

void handle_input(int argc, char *argv[]);

bool is_in_exception_list(const int config_nr);

void SU3_load_config(const unsigned config_nr, mtrx_3x3 *U),
     SU3_write_config(const unsigned config_nr, mtrx_3x3 *U);

void copy_3x3_config(mtrx_3x3 *U, mtrx_3x3 *U_copy);

void SU3_convert_config_in_work(in_cfg_data_type *U_in, work_cfg_data_type *U_work),
     SU3_convert_config_work_out(work_cfg_data_type *U_work, out_cfg_data_type *U_out);

void check_det_1(mtrx_3x3 *U);

void SU3_reunitarize(mtrx_3x3 *U);

/*============================JONIVAR'S CODE===============================*/

void block_swap(int *buffer, size_t length);

void block_swap_double(double *buffer, size_t length);

int byte_swap(void* strip, size_t size, size_t length);

#endif
