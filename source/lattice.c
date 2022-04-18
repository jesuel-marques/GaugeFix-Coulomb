#include <complex.h>
#include <stdio.h>  //	Standard C header files
#include <stdlib.h>
#include <string.h>

#include "../SU3_parameters.h"              //	Simulation parameters
#include "../SU3_gaugefixing_parameters.h"  //	Gauge-fixing specific parameters

#include "SU2_ops.h"  //	SU(2) operations
#include "SU3_ops.h"  //	SU(3) operations

#include "lattice.h"  //	Initialization functions and calculations of
                      //	positions and links on the lattice.

extern char configs_dir_name[];

pos_vec assign_position(const pos_index x, const pos_index y, const pos_index z, const pos_index t){
    //	assigns x, y, z and t to a position vector
    pos_vec position;

    position.i = x;
    position.j = y; 
    position.k = z; 
    position.t = t;

    return position;
}

void print_pos_vec(const pos_vec pos) {
    //	prints a position to the screen

    printf("% d %d %d %d\n", pos.i, pos.j, pos.k, pos.t);
}

pos_vec hop_position_positive(const pos_vec u, const lorentz_index mu) {
    //	Calculates the position immediately forward
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    pos_vec u_plus_muhat;

    switch (mu) {
        case x_index:
            u_plus_muhat.i = ((u.i + 1) % Nxyz);
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case y_index:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = ((u.j + 1) % Nxyz);
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case z_index:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = ((u.k + 1) % Nxyz);
            u_plus_muhat.t = u.t;
            break;

        case t_index:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = ((u.t + 1) % Nt);
            break;

        default:
            puts("mu outside of range\n");
            exit(EXIT_FAILURE);
    }

    return u_plus_muhat;
}

pos_vec hop_position_negative(const pos_vec u, const lorentz_index mu) {
    //	Calculates the position immediately behind
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    pos_vec u_minus_muhat;

    switch (mu) {

        case x_index:
            u_minus_muhat.i = (((u.i - 1) % Nxyz + Nxyz) % Nxyz);
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case y_index:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = (((u.j - 1) % Nxyz + Nxyz) % Nxyz);
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case z_index:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = (((u.k - 1) % Nxyz + Nxyz) % Nxyz);
            u_minus_muhat.t = u.t;            
            break;
        
        case t_index:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = (((u.t - 1) % Nt + Nt) % Nt);
            break;

        default:
            puts("mu outside of range\n");
            exit(EXIT_FAILURE);
    }

    return u_minus_muhat;
}


unsigned short position_is_odd(const pos_vec position) {
    //	If a position is odd returns 1, if even returns 0

    return ((position.i ^ position.j ^ position.k ^ position.t) & 1);

    // if position is odd, then the XOR of the first bit of each element
    // of position must be 1. Take AND with 1 select this first bit.
}

unsigned short position_is_even(const pos_vec position) {
    //	If a position is even returns 1, if odd returns 0

    return !((position.i ^ position.j ^ position.k ^ position.t) & 1);

    // if position is even, then the XOR of the first bit of each element
    // of position must be 0. Taking AND with 1 select this first bit. Take the NOT of 
    // the odd code, because want 1 for even and 0 for odd.
}


void test_allocation(const void * pointer, const char * location ){ 
    //	Test if allocation was successful.
    if ( pointer == NULL ) {
			
			printf("Memory allocation failed at %s.\n", location);
			exit(EXIT_FAILURE); 
	}
}

matrix_3x3_double *get_link(matrix_3x3_double *U, const pos_vec position, const lorentz_index mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + (((((position.t * Nxyz + position.k) * Nxyz + position.j) * Nxyz + position.i) * d + mu));
}

matrix_3x3_float *get_link_f(matrix_3x3_float *U, const pos_vec position, const lorentz_index mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + (((((position.t * Nxyz + position.k) * Nxyz + position.j) * Nxyz + position.i) * d + mu));
}

void get_link_matrix(matrix_3x3_double * U, const pos_vec position, const lorentz_index mu, direction dir, matrix_3x3_double * u){
	
	// Gets forward or backward link at given position and mu
	// and copies it to u.

    matrix_3x3_double u_aux;

	if (dir == FRONT) {

		SU3_copy(get_link(U, position, mu), u);
		//	Link in the positive way is what is stored in U

	}
	else if (dir == REAR) {

		SU3_hermitean_conjugate(get_link(U, hop_position_negative(position, mu), mu), u);
		//	U_(-mu)(n)=(U_mu(n-mu))^\dagger

	}
}

char *name_configuration_file(const unsigned config) {
    char configs_dir_name_local[max_length_name];
    strcpy(configs_dir_name_local, configs_dir_name);
    
    char config_filename[max_length_name];
    //sprintf(config_filename, "NewFormConfig_%d_beta_5.700_Nxyz_%d_Nt_%d.txt", config, Nxyz, Nt);
    sprintf(config_filename,"Gen2_24x16_%d.cfg",config);
    
    return strcat(configs_dir_name_local, config_filename);
}

void SU3_load_config(const char filename[max_length_name], in_cfg_data_type *U) {
    //	Loads a link configuration from the file with filename to U.

    FILE *config_file;

    // printf("Loading: %s.\t", filename);

    config_file = fopen(filename, "r");

    if (fread(U, sizeof(in_cfg_data_type),  Volume * d , config_file) == Volume * d) {
        puts("U Loaded OK\n");
    } else {
        puts(" Configuration loading failed.\n");
        exit(1);
    }

    fclose(config_file);
}

void SU3_print_config(char filename[max_length_name], const char modifier[max_length_name], out_cfg_data_type *U) {
    //  Loads a link configuration from the file with filename to U.

    FILE *config_file;

    printf("Creating: %s.\t", strcat(filename, modifier));

    config_file = fopen(filename, "w+");

    if (fwrite(U, Nc * Nc * sizeof(out_cfg_data_type), Volume * d , config_file) == Volume * d) {
        // puts("U written OK.\n");
    } else {
        puts(" Configuration writing failed.\n");
    }

    fclose(config_file);
}

void SU3_copy_config(matrix_3x3_double *U, matrix_3x3_double *U_copy) {
    // Copies configuration with pointer U to the one with pointer U_copy.

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < Nxyz; position.k++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.i = 0; position.i < Nxyz; position.i++) {
                        for (lorentz_index mu = 0; mu < d; mu++) {

                            SU3_copy(get_link(U, position, mu), get_link(U_copy, position, mu));
                        
                        }
                    }
                }
            }
        }
}

void SU3_convert_config_fd(matrix_3x3_float *U_float, matrix_3x3_double *U_double) {

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < Nxyz; position.k++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.i = 0; position.i < Nxyz; position.i++) {
                        for (lorentz_index mu = 0; mu < d; mu++) {

                            SU3_convert_fd(get_link_f(U_float, position, mu), get_link(U_double, position, mu));
                        
                        }
                    }
                }
            }
        }
}

void SU3_convert_config_df(matrix_3x3_double *U_double, matrix_3x3_float *U_float) {

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < Nxyz; position.k++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.i = 0; position.i < Nxyz; position.i++) {
                        for (lorentz_index mu = 0; mu < d; mu++) {

                            SU3_convert_df(get_link(U_double, position, mu), get_link_f(U_float, position, mu));
                        
                        }
                    }
                }
            }
        }
}


void SU3_reunitarize(matrix_3x3_double *U) {
    // Reunitarizes the configuration

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < Nxyz; position.k++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.i = 0; position.i < Nxyz; position.i++) {
                        for (lorentz_index mu = 0; mu < d; mu++) {

                            SU3_projection(get_link(U, position, mu));

                        }
                    }
                }
            }
        }
}


/*============================JONIVAR'S CODE===============================*/

void block_swap(int *buffer, size_t length)
{
  size_t i;
  union swapper {
    int integer;
    char pos[4];
  } a, b;

  for( i=0 ; i < length ; i++) {
    a.integer = *buffer;
    b.pos[0] = a.pos[3];
    b.pos[1] = a.pos[2];
    b.pos[2] = a.pos[1];
    b.pos[3] = a.pos[0];
    *buffer  = b.integer;
    buffer++;
  }
}

void block_swap_double(double *buffer, size_t length)
{
  size_t i;
  union swapper {
    double double_number;
    char pos[8];
  } a, b;

  for( i=0 ; i < length ; i++) {
    a.double_number = *buffer;
    b.pos[0] = a.pos[7];
    b.pos[1] = a.pos[6];
    b.pos[2] = a.pos[5];
    b.pos[3] = a.pos[4];
    b.pos[4] = a.pos[3];
    b.pos[5] = a.pos[2];
    b.pos[6] = a.pos[1];
    b.pos[7] = a.pos[0];
    *buffer  = b.double_number;
    buffer++;
  }
}



int byte_swap(void* strip, size_t size, size_t length)
{
  switch (size) {
  case sizeof(float):   /* = 4 */
    block_swap(strip, length / size);
    break;
  case sizeof(double):  /* = 8 */
    block_swap_double(strip, length / size);
    break;
  case 1:
    break;
  default:  /* ERROR: unknown size!! */
    return -1;
  }

  return 0;
}