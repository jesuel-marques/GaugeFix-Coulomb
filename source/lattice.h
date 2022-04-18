#include "../SU3_parameters.h"

#ifndef LATTICE_H
#define LATTICE_H

typedef unsigned short pos_index;

typedef struct {
    pos_index i, j, k;
    pos_index t;
} pos_vec;  //	struct for position vectors


typedef unsigned short lorentz_index;

typedef enum {REAR, FRONT} direction; 

#define x_index 0
#define y_index 1
#define z_index 2
#define t_index 3

typedef float complex in_cfg_data_type;
typedef float complex out_cfg_data_type;


pos_vec assign_position(const pos_index x, const pos_index y, const pos_index z, const pos_index t);

void print_pos_vec(const pos_vec u);

inline pos_vec hop_position_positive(const pos_vec u, const lorentz_index mu);

inline pos_vec hop_position_negative(const pos_vec u, const lorentz_index mu);

inline unsigned short position_is_even(const pos_vec position);

inline unsigned short position_is_odd(const pos_vec position);

void test_allocation(const void * pointer, const char * location );

inline double complex * get_link(double complex *U, const pos_vec position, const lorentz_index mu);

inline float complex *get_link_f(float complex *U, const pos_vec position, const lorentz_index mu);

void get_link_matrix(double complex * U, const pos_vec position, const lorentz_index mu, direction dir, double complex * u);


char * name_configuration_file(const unsigned config);

void SU3_load_config(const char filename[max_length_name], in_cfg_data_type *U);

void SU3_print_config(char filename[max_length_name], const char modifier[max_length_name], out_cfg_data_type *U);

void SU3_copy_config(double complex *U, double complex *U_copy);

void SU3_convert_config_fd(float complex *U_float, double complex *U_double);

void SU3_convert_config_df(double complex *U_double, float complex *U_float);

void SU3_reunitarize(double complex *U);

/*============================JONIVAR'S CODE===============================*/

void block_swap(int *buffer, size_t length);

void block_swap_double(double *buffer, size_t length);

int byte_swap(void* strip, size_t size, size_t length);

#endif
