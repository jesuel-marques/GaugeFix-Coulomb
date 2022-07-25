
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



inline bool position_valid(pos_vec position){
    if (position.t >= 0 && position.t < N_T   &&
        position.i >= 0 && position.i < N_SPC &&
        position.j >= 0 && position.j < N_SPC &&
        position.k >= 0 && position.k < N_SPC    ){
        return true;
    }
    else{
        printf("Position invalid: ");
        print_pos_vec(position);
        printf("\n");
        return false;
    }
}

inline bool position_mu_valid(pos_vec position, lorentz_idx mu){
    if (position.t >= 0 && position.t < N_T   &&
        position.i >= 0 && position.i < N_SPC &&
        position.j >= 0 && position.j < N_SPC &&
        position.k >= 0 && position.k < N_SPC &&
       (mu == T_INDX || mu == X_INDX || 
        mu == Y_INDX || mu == Z_INDX )){

        return true;

    }

    return false;

}

pos_vec assign_position(const pos_index x, const pos_index y, const pos_index z, const pos_index t) {
    //	assigns x, y, z and t to a position vector
    pos_vec position;

    position.i = x;
    position.j = y;
    position.k = z;
    position.t = t;

    if(!position_valid(position)){
        fprintf(stderr, "Position {%d, %d, %d, %d} is invalid. Returning origin instead\n.", x, y, z, t);
        return assign_position(0, 0, 0, 0);
    }

    return position;
}

void print_pos_vec(const pos_vec pos) {
    //	prints a position to the screen

    printf("x: %hu y: %hu z: %hu t: %hu\n", pos.i, pos.j, pos.k, pos.t);
}

inline pos_vec hop_position_positive(const pos_vec u, const lorentz_idx mu) {
    //	Calculates the position immediately forward
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    #ifdef CHECK_POSITION_BOUNDS
        if(!position_mu_valid(u, mu)){
            printf("Position: ");
            print_pos_vec(u);
            printf("\n");
            printf("or mu: %d invalid.\n", mu);
            exit(EXIT_FAILURE);
        }
    #endif

    pos_vec u_plus_muhat;

    unsigned short v;

    switch (mu) {
        case X_INDX:
            u_plus_muhat.i = ((v = u.i + 1) != N_SPC ? v : 0);
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case Y_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = ((v = u.j + 1) != N_SPC ? v : 0);
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = u.t;
            break;

        case Z_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = ((v = u.k + 1) != N_SPC ? v : 0);
            u_plus_muhat.t = u.t;
            break;

        case T_INDX:
            u_plus_muhat.i = u.i;
            u_plus_muhat.j = u.j;
            u_plus_muhat.k = u.k;
            u_plus_muhat.t = ((v = u.t + 1) != N_T   ? v : 0);
            break;

    }

    return u_plus_muhat;
}

inline pos_vec hop_position_negative(const pos_vec u, const lorentz_idx mu) {
    //	Calculates the position immediately behind
    //	in the direction mu, taken into account the
    //	periodic boundary conditions.

    #ifdef CHECK_POSITION_BOUNDS
        if(!position_mu_valid(u, mu)){
            printf("Position: ");
            print_pos_vec(u);
            printf("\n");
            printf("or mu: %d invalid.\n", mu);
            exit(EXIT_FAILURE);
        }
    #endif

    pos_vec u_minus_muhat;

    short v ;

     switch (mu) {
        case X_INDX:
            u_minus_muhat.i = ((v = u.i - 1) != -1 ? v : N_SPC - 1);
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case Y_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = ((v = u.j - 1) != -1 ? v : N_SPC - 1);
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = u.t;
            break;

        case Z_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = ((v = u.k - 1) != -1 ? v : N_SPC - 1);   
            u_minus_muhat.t = u.t;
            break;

        case T_INDX:
            u_minus_muhat.i = u.i;
            u_minus_muhat.j = u.j;
            u_minus_muhat.k = u.k;
            u_minus_muhat.t = ((v = u.t - 1) != -1 ?   v : N_T - 1);
            break;
    }

    return u_minus_muhat;
}

short test_allocation_function(const void *pointer, const char *location) {
    //	Test if allocation was successful.
    if (pointer == NULL) {
        fprintf(stderr, "Memory allocation failed at %s.\n", location);
        return -1;
    }

    return 0;
}


mtrx_3x3 *get_link(mtrx_3x3 *restrict U, const pos_vec position, const lorentz_idx mu) {
    //	Does the pointer arithmetic to get a pointer to link at given position and mu
    #ifdef CHECK_POSITION_BOUNDS
        if(position_mu_valid(position, mu))    
    #endif
            return U + GET_LINK_U(position, mu);
    #ifdef CHECK_POSITION_BOUNDS
        else{
            printf("Program reading config outside of allowed range.\n");
            exit(EXIT_FAILURE);
        }
    #endif
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

    if(position_mu_valid(position, mu)){
        if (dir == FRONT) {
        copy_3x3(get_link(U, position, mu), u);
        //	Link in the positive way is what is stored in U

        } else if (dir == REAR) {
        herm_conj_3x3(get_link(U, hop_position_negative(position, mu), mu), u);
        //	U_(-mu)(n)=(U_mu(n-mu))^\dagger
        }
        else{
            printf("direction is not valid.");
            exit(EXIT_FAILURE); 
        }
    }
    else{
        printf("Program reading config outside of allowed range.\n");
        exit(EXIT_FAILURE);
    }
    
}

mtrx_3x3 *get_gaugetransf(mtrx_3x3 *restrict G, const pos_vec position) {
    //	Does the pointer arithmetic to get a pointer to a gaugea transformation at given position

    #ifdef CHECK_POSITION_BOUNDS
        if(position_valid(position)){
    #endif
            return G + GET_GT(position);
    #ifdef CHECK_POSITION_BOUNDS
        }
        else{
        
            printf("Program reading gauge-transformation outside of allowed position range.\n");
            exit(EXIT_FAILURE);
        }
    #endif
}

in_gt_data_type  *get_gaugetransf_in (in_gt_data_type  * restrict G_in,  const pos_vec position) {
    //	Does the pointer arithmetic to get a pointer to a gaugea transformation at given position
    return G_in + GET_GT(position);
}

out_gt_data_type *get_gaugetransf_out(out_gt_data_type * restrict G_out, const pos_vec position) {
    //	Does the pointer arithmetic to get a pointer to a gaugea transformation at given position
    return G_out + GET_GT(position);
}

// void copy_3x3_config(mtrx_3x3 *U, mtrx_3x3 *U_copy) {
//     // Copies configuration with pointer U to the one with pointer U_copy.

// omp_parallel_for
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

#ifdef CONV_TO_WORKING_PRECISION
void SU3_convert_config_in_work(in_cfg_data_type * restrict U_in, work_mtrx_data_type * restrict U_work) {
    omp_parallel_for
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {
                            convert_in_work_cfg_3x3(get_link_in(U_in, position, mu), get_link(U_work, position, mu));
                        }
                    }
                }
            }
        }
}

#endif

#ifdef CONV_GT_TO_WORKING_PRECISION

 void SU3_convert_gt_in_work(in_gt_data_type * restrict G_in, work_mtrx_data_type * restrict G_work) {
        omp_parallel_for
            // Paralelizing by slicing the time extent
            for (pos_index t = 0; t < N_T; t++) {
                pos_vec position;
                position.t = t;
                for (position.k = 0; position.k < N_SPC; position.k++) {
                    for (position.j = 0; position.j < N_SPC; position.j++) {
                        for (position.i = 0; position.i < N_SPC; position.i++) {                    
                                convert_in_work_gt_3x3(get_gaugetransf_in(G_in, position), get_gaugetransf(G_work, position));                    
                        }
                    }
                }
            }
    }
#endif


#ifdef CONV_FROM_WORKING_PRECISION

void SU3_convert_config_work_out(work_mtrx_data_type * restrict U_work, out_cfg_data_type * restrict U_out) {
    omp_parallel_for
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {
                            convert_work_out_cfg_3x3(get_link(U_work, position, mu), get_link_out(U_out, position, mu));
                        }
                    }
                }
            }
        }
}

#endif

#ifdef CONV_GT_FROM_WORKING_PRECISION

void SU3_convert_gt_work_out(work_mtrx_data_type *G_work, out_gt_data_type *G_out) {
    omp_parallel_for
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {                    
                            convert_work_out_gt_3x3(get_gaugetransf(G_work, position), get_gaugetransf_out(G_out, position));                    
                    }
                }
            }
        }
}

#endif

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

//     omp_parallel_for
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

// void SU3_reunitarize_U_G(mtrx_3x3 *U) {
//     // Reunitarizes the configuration
//     // printf("here\n");
//     loop_over_links(1, projection_SU3, U, DOUBLE);
// }

double average_det(mtrx_3x3 *restrict U) {
    double det = 0.0;
    // Reunitarizes the configuration

    #pragma omp parallel for num_threads(NUM_THREADS) \
                        schedule(dynamic) reduction(+ : det)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {
                            det += determinant_3x3(get_link(U, position, mu));
                        }
                    }
                }
            }
        }

    det /= (DIM * VOLUME);

    return det;
}

short SU3_reunitarize_U_G(mtrx_3x3 *restrict U, mtrx_3x3 *restrict G) {
    // Reunitarizes the configuration
    short exit_status = 0;

    #pragma omp parallel for num_threads(NUM_THREADS) reduction (|:exit_status) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        exit_status |= projection_SU3(get_gaugetransf(G, position));
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {
                            exit_status |= projection_SU3(get_link(U, position, mu));
                        }
                    }
                }
            }
        }

    if(!exit_status)
        fprintf(stderr, "At least one matrix could not be projected to SU(3)\n");

    return exit_status;
}


// short SU3_reunitarize_U(mtrx_3x3 *restrict U) {
//     // Reunitarizes the configuration
//     short exit_status = 0;

//     #pragma omp parallel for num_threads(NUM_THREADS) reduction (|:exit_status) schedule(dynamic)
//         // Paralelizing by slicing the time extent
//         for (pos_index t = 0; t < N_T; t++) {
//             pos_vec position;
//             position.t = t;
//             for (position.k = 0; position.k < N_SPC; position.k++) {
//                 for (position.j = 0; position.j < N_SPC; position.j++) {
//                     for (position.i = 0; position.i < N_SPC; position.i++) {
//                         for (lorentz_idx mu = 0; mu < DIM; mu++) {
//                             exit_status |= projection_SU3(get_link(U, position, mu));
//                         }
//                     }
//                 }
//             }
//         }

//     return exit_status;
// }


// short SU3_reunitarize_G(mtrx_3x3 *restrict G) {
//     // Reunitarizes the configuration
//     short exit_status = 0;

//     #pragma omp parallel for num_threads(NUM_THREADS) reduction (|:exit_status) schedule(dynamic)
//         // Paralelizing by slicing the time extent
//         for (pos_index t = 0; t < N_T; t++) {
//             pos_vec position;
//             position.t = t;
//             for (position.k = 0; position.k < N_SPC; position.k++) {
//                 for (position.j = 0; position.j < N_SPC; position.j++) {
//                     for (position.i = 0; position.i < N_SPC; position.i++) {
//                         exit_status |= projection_SU3(get_gaugetransf(G, position));
//                     }
//                 }
//             }
//         }

//     return exit_status;
// }