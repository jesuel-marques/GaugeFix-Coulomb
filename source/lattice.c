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

void print_pos_vec(const pos_vec pos) {
    //	prints a position to the screen

    printf("% d %d %d %d\n", pos.t, pos.i, pos.j, pos.k);
}

pos_vec hop_position_positive(const pos_vec u, const unsigned short mu) {
    //	Calculates the position immediately forward
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    pos_vec u_plus_muhat;

    switch (mu) {
        case 0:
            u_plus_muhat.t = ((u.t + 1) % Nt);
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            break;

        case 1:
            u_plus_muhat.t = u.t;
            u_plus_muhat.i = ((u.i + 1) % Nxyz);
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            break;

        case 2:
            u_plus_muhat.t = u.t;
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = ((u.j + 1) % Nxyz);
            u_plus_muhat.k = u.k;
            break;

        case 3:
            u_plus_muhat.t = u.t;
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = ((u.k + 1) % Nxyz);
            break;

        default:
            printf("mu outside of range");
            exit(EXIT_FAILURE);
    }

    return u_plus_muhat;
}

pos_vec hop_position_negative(const pos_vec u, const unsigned short mu) {
    //	Calculates the position immediately behind
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    pos_vec u_minus_muhat;

    switch (mu) {
        case 0:
            u_minus_muhat.t = (((u.t - 1) % Nt + Nt) % Nt);
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            break;

        case 1:
            u_minus_muhat.t = u.t;
            u_minus_muhat.i = (((u.i - 1) % Nxyz + Nxyz) % Nxyz);
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            break;

        case 2:
            u_minus_muhat.t = u.t;
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = (((u.j - 1) % Nxyz + Nxyz) % Nxyz);
            u_minus_muhat.k = u.k;
            break;

        case 3:
            u_minus_muhat.t = u.t;
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = (((u.k - 1) % Nxyz + Nxyz) % Nxyz);
            break;

        default:
            printf("mu outside of range");
            exit(EXIT_FAILURE);
    }

    return u_minus_muhat;
}


unsigned short position_is_odd(const pos_vec position) {
    //	If a position is odd returns 1, if even returns 0

    return ((position.t ^ position.i ^ position.j ^ position.k) & 1);

    // if position is odd, then the XOR of the first bit of each element
    // of position must be 1. Take AND with 1 select this first bit.
}

unsigned short position_is_even(const pos_vec position) {
    //	If a position is even returns 1, if odd returns 0

    return !((position.t ^ position.i ^ position.j ^ position.k) & 1);

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

double complex *get_link(double complex *U, const pos_vec position, const unsigned short mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + (((((position.t * Nxyz + position.i) * Nxyz + position.j) * Nxyz + position.k) * d + mu) * 3 * 3);
}

float complex *get_link_f(float complex *U, const pos_vec position, const unsigned short mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + (((((position.t * Nxyz + position.i) * Nxyz + position.j) * Nxyz + position.k) * d + mu) * 3 * 3);
}

void get_link_matrix(double complex * U, pos_vec position, int mu, int direction, double complex * u){
	
	// Gets forward or backward link at given position and mu
	// and copies it to u.

	if (direction == 1) {

		SU3_copy(get_link(U, position, mu), u);
		//	Link in the positive way is what is stored in U

	}
	else if (direction == -1) {

		SU3_hermitean_conjugate(get_link(U, hop_position_negative(position, mu), mu), u);
		//	U_(-mu)(n)=(U_mu(n-mu))^\dagger

	}
}

char *name_configuration_file(const unsigned config) {
    char configs_dir_name_local[max_length_name];
    strcpy(configs_dir_name_local, configs_dir_name);
    
    char config_filename[max_length_name];
    //sprintf(config_filename, "NewFormConfig_%d_beta_5.700_Nxyz_%d_Nt_%d.txt", config, Nxyz, Nt);
    sprintf(config_filename,"lend_Gen2_24x16_%d.cfg",config);
    
    return strcat(configs_dir_name_local, config_filename);
}

void SU3_load_config(const char filename[max_length_name], float complex *U) {
    //	Loads a link configuration from the file with filename to U.

    FILE *config_file;

    //printf("Loading: %s.\t", filename);

    config_file = fopen(filename, "r");

    if (fread(U, 3 * 3 * sizeof(float complex),  Volume * d , config_file)) {
        //printf("U Loaded OK\n");
    } else {
        printf(" Configuration loading failed.\n");
        exit(1);
    }

    fclose(config_file);
}

void SU3_print_config(char filename[max_length_name], const char modifier[max_length_name], double complex *U) {
    //  Loads a link configuration from the file with filename to U.

    FILE *config_file;

    printf("Creating: %s.\t", strcat(filename, modifier));

    config_file = fopen(filename, "w+");

    if (fwrite(U, Volume * d * 3 * 3 * sizeof(double complex), 1, config_file) == 1) {
        printf("U written OK.\n");
    } else {
        printf(" Configuration writing failed.\n");
    }

    fclose(config_file);
}

void SU3_print_config_f(char filename[max_length_name], const char modifier[max_length_name], float complex *U) {
    //  Loads a link configuration from the file with filename to U.

    FILE *config_file;

    printf("Creating: %s.\t", strcat(filename, modifier));

    config_file = fopen(filename, "w+");

    if (fwrite(U, Volume * d * 3 * 3 * sizeof(float complex), 1, config_file) == 1) {
        printf("U written OK.\n");
    } else {
        printf(" Configuration writing failed.\n");
    }

    fclose(config_file);
}

void SU3_copy_config(double complex *U, double complex *U_copy) {
    // Copies configuration with pointer U to the one with pointer U_copy.

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (unsigned short t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.i = 0; position.i < Nxyz; position.i++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.k = 0; position.k < Nxyz; position.k++) {
                        for (unsigned short mu = 0; mu < d; mu++) {

                            SU3_copy(get_link(U, position, mu), get_link(U_copy, position, mu));
                        
                        }
                    }
                }
            }
        }
}

void SU3_convert_config_fd(float complex *U_float, double complex *U_double) {

    double complex u[3][3];

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (unsigned short t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.i = 0; position.i < Nxyz; position.i++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.k = 0; position.k < Nxyz; position.k++) {
                        for (unsigned short mu = 0; mu < d; mu++) {

                            SU3_convert_fd(get_link_f(U_float, position, mu), get_link(U_double, position, mu));
                        
                        }
                    }
                }
            }
        }
}

void SU3_convert_config_df(double complex *U_double, float complex *U_float) {

    double complex u[3][3];

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (unsigned short t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.i = 0; position.i < Nxyz; position.i++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.k = 0; position.k < Nxyz; position.k++) {
                        for (unsigned short mu = 0; mu < d; mu++) {

                            SU3_convert_df(get_link(U_double, position, mu), get_link_f(U_float, position, mu));
                        
                        }
                    }
                }
            }
        }
}


void SU3_reunitarize(double complex *U) {
    // Reunitarizes the configuration

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (unsigned short t = 0; t < Nt; t++) {
            pos_vec position;
            position.t = t;
            for (position.i = 0; position.i < Nxyz; position.i++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.k = 0; position.k < Nxyz; position.k++) {
                        for (unsigned short mu = 0; mu < d; mu++) {

                            SU3_projection(get_link(U, position, mu));

                        }
                    }
                }
            }
        }
}
