#include "../SU3_parameters.h"   

#ifndef LATTICE_H
#define LATTICE_H

typedef struct {
    unsigned short t;
    unsigned short i, j, k;
} pos_vec;  //	struct for position vectors

pos_vec add_pos_vec(const pos_vec u, const pos_vec v);

void print_pos_vec(const pos_vec u);

pos_vec hop_position_positive(const pos_vec u, const unsigned short mu);

pos_vec hop_position_negative(const pos_vec u, const unsigned short mu);

unsigned short position_is_even(const pos_vec position);

unsigned short position_is_odd(const pos_vec position);

void test_allocation(const void * pointer, const char * location );

double complex * get_link(double complex *U, const pos_vec position, const unsigned short mu);

float complex *get_link_f(float complex *U, const pos_vec position, const unsigned short mu);

void get_link_matrix(double complex * U, pos_vec position, int mu, int direction, double complex * u);


char * name_configuration_file(const unsigned config);

void SU3_load_config(const char filename[max_length_name], float complex *U);

void SU3_print_config(char filename[max_length_name], const char modifier[max_length_name], double complex *U);

void SU3_print_config_f(char filename[max_length_name], const char modifier[max_length_name], float complex *U);

void SU3_copy_config(double complex *U, double complex *U_copy);

void SU3_convert_config_fd(float complex *U_float, double complex *U_double);

void SU3_convert_config_df(double complex *U_double, float complex *U_float);

void SU3_reunitarize(double complex *U);

#endif
