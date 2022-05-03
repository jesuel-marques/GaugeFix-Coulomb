#include <complex.h>
#include <stdio.h>  //	Standard C header files
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <stdarg.h>

#include "../SU3_parameters.h"              //	Simulation parameters
#include "../SU3_gaugefixing_parameters.h"  //	Gauge-fixing specific parameters

#include "SU2_ops.h"  //	SU(2) operations
#include "SU3_ops.h"  //	SU(3) operations

#include "lattice.h"  //	Initialization functions and calculations of
                      //	positions and links on the lattice.

extern char config_template[];
extern char configs_dir_name[];
extern char extension_in[];
extern char extension_out[];

extern const int config_exception_list[];

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

inline pos_vec hop_position_positive(const pos_vec u, const lorentz_idx mu) {
    //	Calculates the position immediately forward
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    pos_vec u_plus_muhat;

    switch (mu) {
        case X_INDX:
            u_plus_muhat.i = ((u.i + 1) % N_SPC);
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case Y_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = ((u.j + 1) % N_SPC);
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case Z_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = ((u.k + 1) % N_SPC);
            u_plus_muhat.t = u.t;
            break;

        case T_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = ((u.t + 1) % N_T);
            break;

        default:
            fprintf(stderr, "mu outside of range\n");
            exit(EXIT_FAILURE);
    }

    return u_plus_muhat;
}

inline pos_vec hop_position_negative(const pos_vec u, const lorentz_idx mu) {
    //	Calculates the position immediately behind
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    pos_vec u_minus_muhat;

    switch (mu) {

        case X_INDX:
            u_minus_muhat.i = (((u.i - 1) % N_SPC + N_SPC) % N_SPC);
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case Y_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = (((u.j - 1) % N_SPC + N_SPC) % N_SPC);
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case Z_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = (((u.k - 1) % N_SPC + N_SPC) % N_SPC);
            u_minus_muhat.t = u.t;            
            break;
        
        case T_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = (((u.t - 1) % N_T + N_T) % N_T);
            break;

        default:
            fprintf(stderr, "mu outside of range\n");
            exit(EXIT_FAILURE);
    }

    return u_minus_muhat;
}

void test_allocation_in_function(const void * pointer, const char * location ){ 
    //	Test if allocation was successful.
    if ( pointer == NULL ) {
			
			fprintf(stderr, "Memory allocation failed at %s.\n", location);
			exit(EXIT_FAILURE); 
	}
}

mtrx_3x3 *get_link(mtrx_3x3 * restrict U, const pos_vec position, const lorentz_idx mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + GET_LINK(position, mu);
}

in_cfg_data_type *get_link_in(in_cfg_data_type *U, const pos_vec position, const lorentz_idx mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + GET_LINK(position, mu);
}

out_cfg_data_type *get_link_out(out_cfg_data_type *U, const pos_vec position, const lorentz_idx mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + GET_LINK(position, mu);
}

void get_link_matrix(mtrx_3x3 * restrict  U, const pos_vec position, const lorentz_idx mu, direction dir, mtrx_3x3 * restrict  u){
	
	// Gets forward or backward link at given position and mu
	// and copies it to u.

    mtrx_3x3 u_aux;

	if (dir == FRONT) {

		copy_3x3(get_link(U, position, mu), u);
		//	Link in the positive way is what is stored in U

	}
	else if (dir == REAR) {

		SU3_herm_conj(get_link(U, hop_position_negative(position, mu), mu), u);
		//	U_(-mu)(n)=(U_mu(n-mu))^\dagger

	}
}

void greeter_function(char * program_name){

    printf("Program %s compiled at %s on %s\n", program_name, __TIME__, __DATE__);
    printf("Using C version: %ld\n", __STDC_VERSION__);
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

bool is_in_exception_list(const int config_nr){

    for( int i = 0; config_exception_list[i] != -1; i++ ){

        if(config_nr ==  config_exception_list[i]){
            
            return true;
        }
    
    }

    return false;
}

// int extract_config(const unsigned config_nr){   
//     char config_filename[MAX_LENGTH_NAME];

//     char command_lime[MAX_LENGTH_NAME];

//     name_configuration_file(config_nr, config_filename);

//     strcat(config_filename, extension_in);

//     sprintf(command_lime, "lime_extract_record /data/gamow/hot/Gen2/%dx%d/Run1_cfg_%d.lime 2 4 %s", N_SPC, N_T, config_filename);
    
//     if(system(command_lime) != 0){
//         printf("Problem extracting config");
//         exit(1);
//     }

//     return 0;
// }

const char *name_configuration_file(const unsigned config_nr, char * config_filename) {
    
    //sprintf(config_filename, "NewFormConfig_%d_beta_5.700_Nxyz_%d_Nt_%d.txt", config, N_SPC, N_T);
    sprintf(config_filename,"%s%s_%dx%d_%d", configs_dir_name, config_template, N_SPC, N_T, config_nr);
}

void SU3_load_config(const unsigned config_nr, mtrx_3x3 *U) {
    //	Loads a link configuration from the file with filename to U.

    in_cfg_data_type * U_in;

    #ifdef NEED_CONV_TO_WORKING_PRECISION
        
        U_in = (in_cfg_data_type *) malloc(VOLUME * DIM * sizeof(in_cfg_data_type));
	    TEST_ALLOCATION(U_in);

    #else

        U_in = (in_cfg_data_type *) U;

    #endif

    char config_filename[MAX_LENGTH_NAME];
	name_configuration_file(config_nr, config_filename);

    FILE *config_file;

    strcat(config_filename, extension_in);

    printf("Loading: %s.\n", config_filename);
    if((config_file = fopen(config_filename, "rb")) == NULL){
        
        fprintf(stderr, "Error opening config from file %s.\n", config_filename);
        exit(EXIT_FAILURE);

    }


    if (fread(U_in, sizeof(in_cfg_data_type),  VOLUME * DIM , config_file) == VOLUME * DIM) {
        printf("U Loaded OK for config %d.\n", config_nr);
    } else {
        fprintf(stderr, "Configuration loading failed for config %d.\n", config_nr);
        exit(EXIT_FAILURE);
    }

    // fseek(config_file, 1, SEEK_CUR);f

    fgetc(config_file);

    if( !feof(config_file) ) {
        fprintf(stderr,"File has not been read till the end. Check lattice sizes.\n");
        exit(EXIT_FAILURE);    
    }


    fclose(config_file);

    

    #ifdef NEED_BYTE_SWAP_IN

    	byte_swap(U_in, sizeof(in_data_type) / 2 , VOLUME * DIM * sizeof(in_cfg_data_type));
    
    #endif
    
    #ifdef NEED_CONV_TO_WORKING_PRECISION
        SU3_convert_config_in_work(U_in, U);    
	    free(U_in);
    #endif
}

void SU3_write_config(const unsigned config_nr, mtrx_3x3 *U) {
    //  Loads a link configuration from the file with filename to U.
    
    out_cfg_data_type * U_out;

    #ifdef NEED_CONV_FROM_WORKING_PRECISION
        
        U_out = (out_cfg_data_type *) malloc(VOLUME * DIM * sizeof(out_cfg_data_type));
	    TEST_ALLOCATION(U_out);

     	SU3_convert_config_work_out(U, U_out);

    #else

        U_out = (out_cfg_data_type *) U;

    #endif

    #ifdef NEED_BYTE_SWAP_OUT

    	byte_swap(U_out, sizeof(out_data_type) / 2 , VOLUME * DIM * sizeof(out_cfg_data_type));
    
    #endif

    char config_filename[MAX_LENGTH_NAME];
	name_configuration_file(config_nr, config_filename);
    
    FILE *config_file;

    printf("Creating: %s.\n", strcat(config_filename, extension_out));
    if((config_file = fopen(config_filename, "wb")) == NULL){
        
        fprintf(stderr, "Error creating file %s for config.\n", config_filename);
        exit(EXIT_FAILURE);

    }

    if (fwrite(U_out, sizeof(out_cfg_data_type), VOLUME * DIM , config_file)
                                                                 == VOLUME * DIM) {
       
        printf("U written OK for config %d.\n", config_nr);

    } else {

        fprintf(stderr,"Configuration writing failed for config %d.\n", config_nr);

    }

    fclose(config_file);

    #ifdef NEED_CONV_FROM_DOUBLE        
        free(U_out);
    #endif
    
}


void copy_3x3_config(mtrx_3x3 *U, mtrx_3x3 *U_copy) {
    // Copies configuration with pointer U to the one with pointer U_copy.

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {

                            copy_3x3(get_link(U, position, mu), get_link(U_copy, position, mu));
                        
                        }
                    }
                }
            }
        }
}

void SU3_convert_config_in_work(in_cfg_data_type *U_in, work_cfg_data_type *U_work) {

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {

                            convert_in_work_3x3(get_link_in(U_in, position, mu), get_link(U_work, position, mu));
                        
                        }
                    }
                }
            }
        }
}

void SU3_convert_config_work_out(work_cfg_data_type *U_work, out_cfg_data_type *U_out) {

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {

                            convert_work_out_3x3(get_link(U_work, position, mu), get_link_out(U_out, position, mu));
                        
                        }
                    }
                }
            }
        }
}

// typedef void (*loop_f)(mtrx_3x3 *, ...);

// typedef enum {
//     CHAR,
//     INT,
//     FLOAT,
//     DOUBLE
// } TYPE;

// void loop_over_links(int num_arg, loop_f f, mtrx_3x3 *U, TYPE type,  ...) {
//     // Generic loop

//     va_list arguments;
//     va_start(arguments, num_arg+1);

//     #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
//         // Paralelizing by slicing the time extent
//         for (pos_index t = 0; t < N_T; t++) {
//             pos_vec position;
//             position.t = t;
//             for (position.k = 0; position.k < N_SPC; position.k++) {
//                 for (position.j = 0; position.j < N_SPC; position.j++) {
//                     for (position.i = 0; position.i < N_SPC; position.i++) {
//                         for (lorentz_idx mu = 0; mu < DIM; mu++) {
//                             switch(type){
//                                 case FLOAT:
//                                     f(get_link(U, position, mu), get_link_f(va_arg(arguments, mtrx_3x3_float *),position, mu));
//                                     break;
//                                 case DOUBLE:
//                                     f(get_link(U, position, mu), get_link(va_arg(arguments, mtrx_3x3 *),position, mu));
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

// void SU3_reunitarize(mtrx_3x3 *U) {
//     // Reunitarizes the configuration
//     // printf("here\n");
//     loop_over_links(1, projection_SU3, U, DOUBLE);
// }

void check_det_1(mtrx_3x3 * restrict U) {

    double det = 0.0;
    // Reunitarizes the configuration

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) reduction (+:det)
        // Paralelizing by slicing the time extent
        for (unsigned short t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (unsigned short mu = 0; mu < DIM; mu++) {
                            det += determinant_3x3(get_link(U, position, mu));
                        }
                    }
                }
            }
        }

    printf("average determinant: %.15lf\n", det / (DIM * VOLUME));
}

void SU3_reunitarize(mtrx_3x3 * restrict U) {
    // Reunitarizes the configuration

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (unsigned short t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (unsigned short mu = 0; mu < DIM; mu++) {
                            // print_matrix_3x3(get_link(U, position, mu),"link");
                            projection_SU3(get_link(U, position, mu));
                            // print_matrix_3x3(get_link(U, position, mu),"link");
                            // getchar();

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