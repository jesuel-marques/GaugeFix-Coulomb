#include <complex.h>
#include <stdio.h>  //	Standard C header files
#include <stdlib.h>
#include <string.h>

#include <stdarg.h>

#include "../SU3_parameters.h"              //	Simulation parameters
#include "../SU3_gaugefixing_parameters.h"  //	Gauge-fixing specific parameters

#include "SU2_ops.h"  //	SU(2) operations
#include "SU3_ops.h"  //	SU(3) operations

#include "lattice.h"  //	Initialization functions and calculations of
                      //	positions and links on the lattice.

extern char config_template[];
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

    printf("%d %d %d %d\n", pos.i, pos.j, pos.k, pos.t);
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
            fprintf(stderr, "mu outside of range\n");
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
            fprintf(stderr, "mu outside of range\n");
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
			
			fprintf(stderr, "Memory allocation failed at %s.\n", location);
			exit(EXIT_FAILURE); 
	}
}

matrix_3x3_double *get_link(matrix_3x3_double *U, const pos_vec position, const lorentz_index mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + (((((position.t * Nxyz + position.k) * Nxyz + position.j) * Nxyz + position.i) * DIM + mu));
}

matrix_3x3_float *get_link_f(matrix_3x3_float *U, const pos_vec position, const lorentz_index mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + (((((position.t * Nxyz + position.k) * Nxyz + position.j) * Nxyz + position.i) * DIM + mu));
}

void get_link_matrix(matrix_3x3_double * U, const pos_vec position, const lorentz_index mu, direction dir, matrix_3x3_double * u){
	
	// Gets forward or backward link at given position and mu
	// and copies it to u.

    matrix_3x3_double u_aux;

	if (dir == FRONT) {

		copy_3x3(get_link(U, position, mu), u);
		//	Link in the positive way is what is stored in U

	}
	else if (dir == REAR) {

		SU3_hermitean_conjugate(get_link(U, hop_position_negative(position, mu), mu), u);
		//	U_(-mu)(n)=(U_mu(n-mu))^\dagger

	}
}

void handle_input(int argc, char *argv[]){
    if( argc != 3){
		fprintf(stderr, "Input names for config directory and template");
        
        if(argc > 3){
	    	fprintf(stderr, " only");
	    }
	}

    if( argc != 3){
		fprintf(stderr, ".\n");
		exit(EXIT_FAILURE);
	}
    

	strcpy(configs_dir_name, argv[1]);
	strcpy(config_template, argv[2]); 
}


const char *name_configuration_file(const unsigned config_number, char * config_filename) {
    
    //sprintf(config_filename, "NewFormConfig_%d_beta_5.700_Nxyz_%d_Nt_%d.txt", config, Nxyz, Nt);
    sprintf(config_filename,"%s%s_%dx%d_%d", configs_dir_name, config_template, Nxyz, Nt, config_number);
}

void SU3_load_config(const unsigned config_nr, matrix_3x3_double *U) {
    //	Loads a link configuration from the file with filename to U.


    char config_filename[MAX_LENGTH_NAME];
	name_configuration_file(config_nr, config_filename);
    // const char * extension_in = "_clmbgf.cfg";
    const char * extension_in = ".cfg";

    FILE *config_file;

    strcat(config_filename, extension_in);

    printf("Loading: %s.\n", config_filename);
    if((config_file = fopen(config_filename, "rb")) == NULL){
        
        fprintf(stderr, "Error loading config from file %s.\n", config_filename);
        exit(EXIT_FAILURE);

    }

    in_cfg_data_type * U_in = (in_cfg_data_type *) malloc(VOLUME * DIM * sizeof(in_cfg_data_type));
	test_allocation(U_in, "SU3_load_config");

    if (fread(U_in, sizeof(in_cfg_data_type),  VOLUME * DIM , config_file) == VOLUME * DIM) {
        puts("U Loaded OK\n");
    } else {
        fprintf(stderr, " Configuration loading failed.\n");
        exit(EXIT_FAILURE);
    }

    fclose(config_file);

    if(NEED_BYTE_SWAP){
    	byte_swap(U_in, sizeof(float) , VOLUME * DIM * sizeof(in_cfg_data_type));
    }

    SU3_convert_config_fd(U_in, U);
	free(U_in);
}

void SU3_write_config(const unsigned config_nr, matrix_3x3_double *U) {
    //  Loads a link configuration from the file with filename to U.
    
    out_cfg_data_type * U_out = (out_cfg_data_type *) malloc(VOLUME * DIM * sizeof(out_cfg_data_type));
	test_allocation(U_out, "SU3_write_config");

	SU3_convert_config_df(U, U_out);

    char config_filename[MAX_LENGTH_NAME];
	name_configuration_file(config_nr, config_filename);
    const char * extension_out = "_clmbgf.cfg";
    
    FILE *config_file;

    printf("Creating: %s.\n", strcat(config_filename, extension_out));
    if((config_file = fopen(config_filename, "wb")) == NULL){
        
        fprintf(stderr, "Error creating file %s for config.\n", config_filename);
        exit(EXIT_FAILURE);

    }

    if (fwrite(U_out, sizeof(out_cfg_data_type), VOLUME * DIM , config_file) == VOLUME * DIM) {
       
        puts("U written OK.\n");

    } else {

        fprintf(stderr,"Configuration writing failed.\n");

    }

    fclose(config_file);
    free(U_out);
}


void copy_3x3_config(matrix_3x3_double *U, matrix_3x3_double *U_copy) {
    // Copies configuration with pointer U to the one with pointer U_copy.

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < Nxyz; position.k++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.i = 0; position.i < Nxyz; position.i++) {
                        for (lorentz_index mu = 0; mu < DIM; mu++) {

                            copy_3x3(get_link(U, position, mu), get_link(U_copy, position, mu));
                        
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
                        for (lorentz_index mu = 0; mu < DIM; mu++) {

                            convert_fd_3x3(get_link_f(U_float, position, mu), get_link(U_double, position, mu));
                        
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
                        for (lorentz_index mu = 0; mu < DIM; mu++) {

                            convert_df_3x3(get_link(U_double, position, mu), get_link_f(U_float, position, mu));
                        
                        }
                    }
                }
            }
        }
}

// typedef void (*loop_f)(matrix_3x3_double *, ...);

// typedef enum {
//     CHAR,
//     INT,
//     FLOAT,
//     DOUBLE
// } TYPE;

// void loop_over_links(int num_arg, loop_f f, matrix_3x3_double *U, TYPE type,  ...) {
//     // Generic loop

//     va_list arguments;
//     va_start(arguments, num_arg+1);

//     #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
//         // Paralelizing by slicing the time extent
//         for (pos_index t = 0; t < Nt; t++) {
//             pos_vec position;
//             position.t = t;
//             for (position.k = 0; position.k < Nxyz; position.k++) {
//                 for (position.j = 0; position.j < Nxyz; position.j++) {
//                     for (position.i = 0; position.i < Nxyz; position.i++) {
//                         for (lorentz_index mu = 0; mu < DIM; mu++) {
//                             switch(type){
//                                 case FLOAT:
//                                     f(get_link(U, position, mu), get_link_f(va_arg(arguments, matrix_3x3_float *),position, mu));
//                                     break;
//                                 case DOUBLE:
//                                     f(get_link(U, position, mu), get_link(va_arg(arguments, matrix_3x3_double *),position, mu));
//                                     break;
//                                 default:
//                                     break;
//                             }

//                         }
//                     }
//                 }
//             }
//         }

//     va_end(arguments);
// }

// void SU3_reunitarize(matrix_3x3_double *U) {
//     // Reunitarizes the configuration
//     // printf("here\n");
//     loop_over_links(1, projection_SU3, U, DOUBLE);
// }

void SU3_reunitarize(matrix_3x3_double *U) {
    // Reunitarizes the configuration

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (unsigned short t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < Nxyz; position.k++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.i = 0; position.i < Nxyz; position.i++) {
                        for (unsigned short mu = 0; mu < DIM; mu++) {

                            projection_SU3(get_link(U, position, mu));

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