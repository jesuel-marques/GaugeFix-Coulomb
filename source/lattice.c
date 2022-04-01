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

pos_vec add_position_vector(const pos_vec u, const  pos_vec v) {
    //	Adds two position vectors v and u, taking into account the periodic boundary conditions,
    //  and puts result in vplusu.

    pos_vec u_plus_v;

    if ((u.t + v.t) >= 0) {  //	time component
        u_plus_v.t = ((u.t + v.t) % Nt);
    } else {
        u_plus_v.t = (((u.t + v.t) % Nt + Nt) % Nt);
    }

    if ((u.i + v.i) >= 0) {  //	x component
        u_plus_v.i = ((u.i + v.i) % Nxyz);
    } else {
        u_plus_v.i = (((u.i + v.i) % Nxyz + Nxyz) % Nxyz);
    }

    if ((u.j + v.j) >= 0) {  //	y component
        u_plus_v.j = ((u.j + v.j) % Nxyz);
    } else {
        u_plus_v.j = (((u.j + v.j) % Nxyz + Nxyz) % Nxyz);
    }

    if ((u.k + v.k) >= 0) {  //	z component
        u_plus_v.k = ((u.k + v.k) % Nxyz);
    } else {
        u_plus_v.k = (((u.k + v.k) % Nxyz + Nxyz) % Nxyz);
    }

    return u_plus_v;
}

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
            exit(1);
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
            exit(1);
    }

    return u_minus_muhat;
}

unsigned short position_is_even(const pos_vec position) {
    //	If a position is even returns 1, if odd returns 0

    return ((position.t + position.i + position.j + position.k + 1) % 2);
}

unsigned short position_is_odd(const pos_vec position) {
    //	If a position is even returns 1, if odd returns 0

    return ((position.t + position.i + position.j + position.k) % 2);
}
void test_allocation(const void * pointer, const char * location ){ 
    //	Test if allocation was successful.
    if ( pointer == NULL ) {
			
			printf("Memory allocation failed at %s.\n", location);
			exit(1); 
	}
}

double complex *get_link(double complex *U, const pos_vec position, const unsigned short mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + (((((position.t * Nxyz + position.i) * Nxyz + position.j) * Nxyz + position.k) * d + mu) * 3 * 3);
}

char *name_configuration_file(const unsigned short config) {
    char configs_dir_name_local[max_length_name];
    strcpy(configs_dir_name_local, configs_dir_name);
    
    char config_filename[max_length_name];
    sprintf(config_filename, "NewFormConfig_%d_beta_5.700_Nxyz_%d_Nt_%d.txt", config, Nxyz, Nt);

    return strcat(configs_dir_name_local, config_filename);
}

void SU3_load_config(const char filename[max_length_name], double complex *U) {
    //	Loads a link configuration from the file with filename to U.

    FILE *config_file;

    printf("Loading: %s\n", filename);

    config_file = fopen(filename, "r");

    if (fread(U, Volume * d * 3 * 3 * sizeof(double complex), 1, config_file) == 1) {
        printf("U Loaded\n");
    } else {
        printf(" Configuration loading failed.\n");
    }

    fclose(config_file);
}

void SU3_copy_config(double complex *U, double complex *U_copy) {
    // Copies configuration with pointer U to the one with pointer U_copy.

    pos_vec position;

    for (position.t = 0; position.t < Nt; position.t++) {
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

void SU3_reunitarize(double complex *U) {
    // Reunitarizes the configuration

    pos_vec position;
    double complex *u = (double complex *)malloc(3 * 3 * sizeof(double complex));
    double complex *sub_x = (double complex *)malloc(3 * 3 * sizeof(double complex));

    for (position.t = 0; position.t < Nt; position.t++) {
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