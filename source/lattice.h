#ifndef LATTICE_H
#define LATTICE_H
#include <stdbool.h>

#include "../SU3_parameters.h"

#define omp_parallel_basic "omp parallel for num_threads(NUM_THREADS) schedule(dynamic)"
#define omp_parallel_for  _Pragma(omp_parallel_basic)
#define omp_parallel_for_reduction(op, x) _Pragma(omp_parallel_basic ## " reduction(" ## op:x ## ")")


typedef unsigned short pos_index;

typedef struct {
    pos_index i, j, k;
    pos_index t;
} pos_vec;  //	struct for position vectors

#define POSITION_IS_ODD(position)    (((position.i) ^ (position.j) ^ (position.k) ^ (position.t)) & 1)
#define POSITION_IS_EVEN(position)   !(POSITION_IS_ODD(position))


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

#define GET_LINK_U(position, mu) ((((position.t * N_SPC \
                                        + position.k) * N_SPC \
                                            + position.j) * N_SPC \
                                                + position.i) * DIM \
                                                    + mu)

//  Does the pointer arithmetic to get the correct index in the gauge transformation

#define GET_GT(position)      (((position.t * N_SPC \
                                        + position.k) * N_SPC \
                                            + position.j) * N_SPC \
                                                + position.i) 
                                                    

//  Data types definitions

typedef complex float  in_data_type; 
typedef complex float  out_data_type;
typedef complex double work_data_type;

typedef struct {
   complex double m[Nc * Nc];
} mtrx_3x3_double; 

typedef struct {
    complex float m[Nc * Nc];
} mtrx_3x3_float;

typedef mtrx_3x3_double work_cfg_data_type;
typedef work_cfg_data_type mtrx_3x3;

#ifdef NEED_CONV_TO_WORKING_PRECISION
typedef mtrx_3x3_float in_cfg_data_type;
#else
typedef work_cfg_data_type in_cfg_data_type;
#endif

#ifdef NEED_CONV_FROM_WORKING_PRECISION
typedef mtrx_3x3_float out_cfg_data_type;
#else
typedef work_cfg_data_type out_cfg_data_type;
#endif



#define TEST_ALLOCATION(a) test_allocation_function(a, __func__) //  used to test if allocation was successful
#define GREETER()   greeter_function(__FILE__)

pos_vec assign_position(const pos_index x, const pos_index y, const pos_index z, const pos_index t);

void print_pos_vec(const pos_vec u);

pos_vec hop_position_positive(const pos_vec u, const lorentz_idx mu),
        hop_position_negative(const pos_vec u, const lorentz_idx mu);


short test_allocation_function(const void * pointer, const char * location );

mtrx_3x3 * get_link(mtrx_3x3 *U, const pos_vec position, const lorentz_idx mu);

in_cfg_data_type  *get_link_in (in_cfg_data_type  * restrict U_in,  const pos_vec position, const lorentz_idx mu);
out_cfg_data_type *get_link_out(out_cfg_data_type * restrict U_out, const pos_vec position, const lorentz_idx mu);

void get_link_matrix(mtrx_3x3 * restrict U, const pos_vec position, const lorentz_idx mu, direction dir, mtrx_3x3 * restrict u);

// void copy_3x3_config(mtrx_3x3 *U, mtrx_3x3 *U_copy);

#ifdef NEED_CONV_TO_WORKING_PRECISION
    void SU3_convert_config_work_out(work_cfg_data_type * restrict U_work, out_cfg_data_type * restrict U_out);
    void SU3_convert_gaugetransf_in_work (in_cfg_data_type * restrict G_in,     work_cfg_data_type * restrict G_work);

#endif

#ifdef NEED_CONV_FROM_WORKING_PRECISION
    void SU3_convert_config_in_work (in_cfg_data_type * restrict U_in,   work_cfg_data_type * restrict U_work);
    void SU3_convert_gaugetransf_work_out(work_cfg_data_type * restrict G_work,  out_cfg_data_type * restrict G_out );
#endif

double check_det_1(mtrx_3x3 * restrict U);

short SU3_reunitarize_U_G(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G);
// short SU3_reunitarize_U(mtrx_3x3 *restrict U) ;
// short SU3_reunitarize_G(mtrx_3x3 *restrict G);


mtrx_3x3 *get_gaugetransf(mtrx_3x3 * restrict G, const pos_vec position);

in_cfg_data_type  *get_gaugetransf_in (in_cfg_data_type  * restrict G_in,  const pos_vec position);
out_cfg_data_type *get_gaugetransf_out(out_cfg_data_type * restrict G_out, const pos_vec position);




#endif
