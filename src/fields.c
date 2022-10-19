#include <stdlib.h>
#include <stdio.h>

#include <fields.h>
#include <types.h>
#include <lattice.h>
#include <SU3_ops.h>

#include <SU3_parameters.h>

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

    #ifdef CHECK_POSITION_BOUNDS
    if(position_mu_valid(position, mu)){
    #endif
        if (dir == FRONT) {
        copy_3x3(get_link(U, position, mu), u);
        //	Link in the positive way is what is stored in U

        } else if (dir == REAR) {
        herm_conj_3x3(get_link(U, hop_pos_minus(position, mu), mu), u);
        //	U_(-mu)(n)=(U_mu(n-mu))^\dagger
        }
        else{
            printf("direction is not valid.");
            exit(EXIT_FAILURE); 
        }
    #ifdef CHECK_POSITION_BOUNDS
    }
    else{
        printf("Program reading config outside of allowed range.\n");
        exit(EXIT_FAILURE);
    }
    #endif
    
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

// OMP_PARALLEL_FOR
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

//     OMP_PARALLEL_FOR
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

    if(exit_status)
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