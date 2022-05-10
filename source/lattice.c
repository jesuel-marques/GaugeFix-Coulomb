
#include <complex.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>  //	Standard C header files
#include <stdlib.h>
#include <string.h>

#include "../SU3_gaugefixing_parameters.h"  //	Gauge-fixing specific parameters
#include "../SU3_parameters.h"              //	Simulation parameters
#include "SU2_ops.h"                        //	SU(2) operations
#include "SU3_ops.h"                        //	SU(3) operations
                                            //	positions and links on the lattice.

#include "lattice.h"  //	Initialization functions and calculations of

pos_vec assign_position(const pos_index x, const pos_index y, const pos_index z, const pos_index t) {
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

void test_allocation_in_function(const void *pointer, const char *location) {
    //	Test if allocation was successful.
    if (pointer == NULL) {
        fprintf(stderr, "Memory allocation failed at %s.\n", location);
        exit(EXIT_FAILURE);
    }
}

mtrx_3x3 *get_link(mtrx_3x3 *restrict U, const pos_vec position, const lorentz_idx mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + GET_LINK_U(position, mu);
}

in_cfg_data_type *get_link_in(in_cfg_data_type *U, const pos_vec position, const lorentz_idx mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + GET_LINK_U(position, mu);
}

out_cfg_data_type *get_link_out(out_cfg_data_type *U, const pos_vec position, const lorentz_idx mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    return U + GET_LINK_U(position, mu);
}

void get_link_matrix(mtrx_3x3 *restrict U, const pos_vec position, const lorentz_idx mu, direction dir, mtrx_3x3 *restrict u) {
    // Gets forward or backward link at given position and mu
    // and copies it to u.

    mtrx_3x3 u_aux;

    if (dir == FRONT) {
        copy_3x3(get_link(U, position, mu), u);
        //	Link in the positive way is what is stored in U

    } else if (dir == REAR) {
        SU3_herm_conj(get_link(U, hop_position_negative(position, mu), mu), u);
        //	U_(-mu)(n)=(U_mu(n-mu))^\dagger
    }
}

mtrx_3x3 *get_gaugetransf(mtrx_3x3 *restrict G, const pos_vec position) {
    //	Does the pointer arithmetic to get a pointer to a gaugea transformation at given position
    return G + GET_GT(position);
}

in_cfg_data_type *get_gaugetransf_in(in_cfg_data_type *restrict G_in, const pos_vec position) {
    //	Does the pointer arithmetic to get a pointer to a gaugea transformation at given position
    return G_in + GET_GT(position);
}

out_cfg_data_type *get_gaugetransf_out(out_cfg_data_type *restrict G_out, const pos_vec position) {
    //	Does the pointer arithmetic to get a pointer to a gaugea transformation at given position
    return G_out + GET_GT(position);
}



// void copy_3x3_config(mtrx_3x3 *U, mtrx_3x3 *U_copy) {
//     // Copies configuration with pointer U to the one with pointer U_copy.

// #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
//     // Paralelizing by slicing the time extent
//     for (pos_index t = 0; t < N_T; t++) {
//         pos_vec position;
//         position.t = t;
//         for (position.k = 0; position.k < N_SPC; position.k++) {
//             for (position.j = 0; position.j < N_SPC; position.j++) {
//                 for (position.i = 0; position.i < N_SPC; position.i++) {
//                     for (lorentz_idx mu = 0; mu < DIM; mu++) {
//                         copy_3x3(get_link(U, position, mu), get_link(U_copy, position, mu));
//                     }
//                 }
//             }
//         }
//     }
// }

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

int check_det_1(mtrx_3x3 *restrict U) {
    double det = 0.0;
    // Reunitarizes the configuration

#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) reduction(+ \
                                                                              : det)
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

    det /= (DIM * VOLUME);

    printf("average determinant: %.15lf\n", det);

    return det;
}

void SU3_reunitarize(mtrx_3x3 *restrict U, mtrx_3x3 *restrict G) {
    // Reunitarizes the configuration

#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
    // Paralelizing by slicing the time extent
    for (unsigned short t = 0; t < N_T; t++) {
        pos_vec position;
        position.t = t;
        for (position.k = 0; position.k < N_SPC; position.k++) {
            for (position.j = 0; position.j < N_SPC; position.j++) {
                for (position.i = 0; position.i < N_SPC; position.i++) {
                    projection_SU3(get_gaugetransf(G, position));
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

void SU3_convert_gaugetransf_in_work(in_cfg_data_type * restrict G_in, work_cfg_data_type * restrict G_work) {
#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
    // Paralelizing by slicing the time extent
    for (pos_index t = 0; t < N_T; t++) {
        pos_vec position;
        position.t = t;
        for (position.k = 0; position.k < N_SPC; position.k++) {
            for (position.j = 0; position.j < N_SPC; position.j++) {
                for (position.i = 0; position.i < N_SPC; position.i++) {                    
                        convert_in_work_3x3(get_gaugetransf_in(G_in, position), get_gaugetransf(G_work, position));                    
                }
            }
        }
    }
}


void SU3_convert_gaugetransf_work_out(work_cfg_data_type *G_work, out_cfg_data_type *G_out) {
#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
    // Paralelizing by slicing the time extent
    for (pos_index t = 0; t < N_T; t++) {
        pos_vec position;
        position.t = t;
        for (position.k = 0; position.k < N_SPC; position.k++) {
            for (position.j = 0; position.j < N_SPC; position.j++) {
                for (position.i = 0; position.i < N_SPC; position.i++) {                    
                        convert_work_out_3x3(get_gaugetransf(G_work, position), get_gaugetransf_out(G_out, position));                    
                }
            }
        }
    }
}
