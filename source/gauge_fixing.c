#include <complex.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>  //	Standard header C files
#include <stdlib.h>

#include "../SU3_gaugefixing_parameters.h"  //	Gauge-fixing specific parameters
#include "../SU3_parameters.h"              //	Simulation parameters
#include "lattice.h"                        //	Initialization functions and
                                            //	calculations of positions and
                                            //	links on the lattice.

#include "SU2_ops.h"           //	SU(2) operations
#include "SU3_ops.h"           //	SU(3) operations
#include "gauge_fixing.h"
#include "fourvector_field.h"  //	Calculation of A_mu(n) and related things
#include "math_ops.h"          //	Math operations

#define sweeps_to_next_measurement(e2)  10 + (unsigned)(INITIAL_SWEEPS_TO_MEASUREMENT_e2 \
                                        * (1.0 - log10((e2)) / log10(TOLERANCE)))

// static void SU3_local_update_U(mtrx_3x3 * restrict U, const pos_vec position, 
//                                                         const mtrx_3x3 * restrict g) {
//     //	Updates U only at a given position

//     mtrx_3x3 g_dagger;

//     for (lorentz_idx mu = 0; mu < DIM; mu++) {
//         //	U'_mu(x)=g(x).U_mu(x).1 for red-black updates

//         accum_left_prod_3x3(g, get_link(U, position, mu));

//         //	U'_mu(x-mu)=1.U_mu(x-mu).g_dagger(x) for red-black updates

//         SU3_herm_conj(g, &g_dagger);

//         accum_right_prod_3x3(get_link(U, hop_position_negative(position, mu), mu), 
//                                                                                &g_dagger);
//     }
// }

static void SU3_local_update_U_G(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G, const pos_vec position, 
                                                        const mtrx_3x3 * restrict g) {
    //	Updates U only at a given position

    accum_left_prod_3x3(g, get_gaugetransf(G, position));

    mtrx_3x3 g_dagger;

    for (lorentz_idx mu = 0; mu < DIM; mu++) {
        //	U'_mu(x)=g(x).U_mu(x).1 for red-black updates

        accum_left_prod_3x3(g, get_link(U, position, mu));

        //	U'_mu(x-mu)=1.U_mu(x-mu).g_dagger(x) for red-black updates

        SU3_herm_conj(g, &g_dagger);

        accum_right_prod_3x3(get_link(U, hop_position_negative(position, mu), mu), 
                                                                               &g_dagger);
    }
}

void SU3_update_global_U_G(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G) {
    //	Updates U on the whole lattice

    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            mtrx_3x3 * g;

            mtrx_3x3 g_dagger_position_plus_mu;

            mtrx_3x3 * u;
            mtrx_3x3 u_updated;
            pos_vec position;

            position.t = t;

            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        g = get_gaugetransf(G,position);
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {
                            //	U'_mu(x)=g(x).U_mu(x).gdagger(x+mu)
                            u = get_link(U, position, mu);

                            SU3_herm_conj(get_gaugetransf(G, hop_position_positive(position, mu)), &g_dagger_position_plus_mu);
                            
                            prod_three_3x3(g, u, &g_dagger_position_plus_mu, &u_updated);
                            
                            copy_3x3(&u_updated, u);

                        }
                    }
                }
            }
        }
}

inline static void SU3_calculate_w(mtrx_3x3 * restrict U, const pos_vec position,
                                                         mtrx_3x3 * restrict w) {
    //	Calculates 	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision, following the notation in hep-lat/9306018
    //	returns result in w.

    set_null_3x3(w);  //	Initializing w(n)=0

    mtrx_3x3 u_dagger_rear;

    // w(n)	calculation

    for (lorentz_idx mu = 0; mu < DIM - 1 ; mu++) {
        //	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision

        accumulate_3x3(get_link(U, position, mu), w);

        SU3_herm_conj(get_link(U, hop_position_negative(position, mu), mu), 
                                                                &u_dagger_rear);

        accumulate_3x3(&u_dagger_rear, w);
    }
}
    
double SU3_calculate_e2(mtrx_3x3 * restrict U) {
    //	Calculates e2 (defined in hep-lat/0301019v2),
    //	used to find out distance to the gauge-fixed situation.
 
    double e2 = 0.0;

    #pragma omp parallel for reduction (+:e2) num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            mtrx_3x3 div_A;
            matrix_SU3_alg div_A_components;
            pos_vec position;

            position.t = t;
            double e2_slice = 0.0;

            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        
                        SU3_divergence_A(U, position, &div_A);
                        decompose_algebra_SU3(&div_A, &div_A_components);

                        for (SU3_alg_idx a = 1; a <= POW2(Nc)-1; a++) {
                            //	Normalized sum of the squares of the color components
                            //  of the divergence of A.

                            e2_slice += POW2(div_A_components.m[a]);
                        }

                    }
                }
            }
            e2 += e2_slice;
        }

    e2 /= (VOLUME);

    return e2;
}

inline static void SU3_update_sub_LosAlamos(mtrx_3x3 * restrict w, submatrix sub) {
    SU3_color_idx a, b;
    #ifdef OLD_VERSION
    mtrx_3x3 sub_update;
    
    set_null_3x3(&sub_update);

    sub_update.m[ELM(2 - sub, 2 - sub)] = 1.0;
    #endif

    a = sub == T ? 1 : 0;
    b = sub == R ? 1 : 2;

    matrix_2x2_ck matrix_SU2;

    matrix_SU2.m[0] =  (creal(w -> m[ELM(a, a)]) 
                      + creal(w -> m[ELM(b, b)]));
    matrix_SU2.m[1] = -(cimag(w -> m[ELM(a, b)]) 
                      + cimag(w -> m[ELM(b, a)]));
    matrix_SU2.m[2] = -(creal(w -> m[ELM(a, b)]) 
                      - creal(w -> m[ELM(b, a)]));
    matrix_SU2.m[3] = -(cimag(w -> m[ELM(a, a)]) 
                      - cimag(w -> m[ELM(b, b)]));

    SU2_projection(&matrix_SU2);
    #ifdef OLD_VERSION
    sub_update.m[ELM(a, a)] =      matrix_SU2.m[0] 
                             + I * matrix_SU2.m[3];
    sub_update.m[ELM(a, b)] =      matrix_SU2.m[2]
                             + I * matrix_SU2.m[1];
    sub_update.m[ELM(b, a)] =     -matrix_SU2.m[2] 
                             + I * matrix_SU2.m[1];
    sub_update.m[ELM(b, b)] =      matrix_SU2.m[0] 
                             - I * matrix_SU2.m[3];

    // print_matrix_3x3(&sub_update,"update su(3)",10);
    // printf("sub: %u a: %u b: %u", sub, a, b);
    // print_matrix_3x3(w, "antes", 10);
    accum_left_prod_3x3(&sub_update, w);
    // print_matrix_3x3(w, "depois", 10);
    #else 
    accum_prod_SU2_3x3(&matrix_SU2, w, sub);
    #endif
}

inline static void SU3_LosAlamos_common_block(mtrx_3x3 * restrict w, 
                                            mtrx_3x3 * restrict total_update) {
    //	Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
    //	algorithm for SU(3), with a division of the update matrix in submatrices
    //	following the Cabbibo-Marinari trick. Actual update is obtained after a number
    //	of "hits" to be performed one after another.

    mtrx_3x3 w_inv_old;

    inverse_3x3(w, &w_inv_old);

    for (unsigned short hits = 1; hits <= MAX_HITS; hits++) {
        //	Each hit contains the Cabbibo-Marinari subdivision
        for (submatrix sub = R; sub <= T; sub++) {
            //	Submatrices are indicated by numbers from 0 to 2

            SU3_update_sub_LosAlamos(w, sub);
            
            //	Updates matrix to total_update. It is the
            //	accumulated updates from the hits.
        }
    }

    prod_3x3(w, &w_inv_old, total_update);
    
}

inline static void SU3_gaugefixing_overrelaxation(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G, const pos_vec position) {
    //	Generalization of the algorithm described in hep-lat/0301019v2, using the
    //	Cabbibo-Marinari submatrices trick.
    //	It updates the g at the given position.

    mtrx_3x3 w;

    SU3_calculate_w(U, position, &w);  //	Calculating w(n)=h(n) for red black subdivision

    mtrx_3x3 update_LA;

    SU3_LosAlamos_common_block(&w, &update_LA);

    /*	The above function determines update_LA which would be the naïve
   	update to bring the local function to its mininum. However,
 	using overrelaxation, which
   	means using update_LA^omega instead of  update_LA, where 1<omega<2 
    makes the converge faster. update_LA^omega is calculated using 
    the first two terms of the binomial expansion: 
    update_LA^omega=I+omega(update_LA-I)+...=I(1-omega)+omega*update_LA+...*/

    mtrx_3x3 update_OR;

    /* update_OR = update_LA^omega = Proj_SU3((I(1-omega)+omega*update_LA) */
    power_3x3_binomial(&update_LA, OMEGA_OR, &update_OR);

    projection_SU3(&update_OR);

    SU3_local_update_U_G(U, G, position, &update_OR);
}

unsigned SU3_gauge_fix(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G, const unsigned short config) {
    //	Fix the gauge and follows the process by calculating e2;

    pos_vec position;

    double e2;

    unsigned sweep = 0;
    unsigned sweeps_to_measurement_e2 = INITIAL_SWEEPS_TO_MEASUREMENT_e2;
    //	Counter to the number of sweeps to fix config to Landau gauge
	e2 = SU3_calculate_e2(U);
	printf("Sweeps in config %5d: %8d. e2: %3.2E \n", config, sweep, e2);
    
    while (1) {

         #pragma omp parallel for num_threads(NUM_THREADS) private(position) schedule(dynamic)
            // Paralelizing by slicing the time extent
            for (pos_index t = 0; t < N_T; t++) {
                position.t = t;
                for (position.k = 0; position.k < N_SPC; position.k++) {
                    for (position.j = 0; position.j < N_SPC; position.j++) {
                        for (position.i = 0; position.i < N_SPC; position.i++) {
                            !((POSITION_IS_EVEN(position) + sweep) % 2) ?
                                //	Implementation of the checkerboard subdivision 
                                //  of the lattice.
                                
                                SU3_gaugefixing_overrelaxation(U, G, position)
                                //  The actual gauge-fixing algorithm
                                                                    
                                : 0;
                            
                        }
                    }
                }
            }

        sweep++;

        !(sweep % (INITIAL_SWEEPS_TO_MEASUREMENT_e2 / 5) ) ?
            printf("Sweeps in config %5d: %8d.\n", config, sweep) : 0; 

        if(sweep == sweeps_to_measurement_e2){
            e2 = SU3_calculate_e2(U);
	        printf("Sweeps in config %5d: %8d. e2: %3.2E \n", config, sweep, e2);            //	Gauge-fixing index, indicates how far we are to the Landau-gauge.
            //  It will be less than the TOLERANCE,
            //	when the gauge-fixing is considered to be attained.
            //	Following the notation of hep-lat/0301019v2
            
            if (e2 <= TOLERANCE) {
                break;

            } else {
                sweeps_to_measurement_e2 += 
                    sweeps_to_next_measurement(e2);
                //	No need to calculate e2 all the time
                //	because it will take some hundreds of sweeps
                //	to fix the gauge.
            }
        }
        

        !(sweep % SWEEPS_TO_REUNITARIZATION) ? SU3_reunitarize(U) : 0;
                                               
    }
    SU3_reunitarize(U);

    printf("Sweeps needed to gauge-fix config %d: %d. e2: %3.2E \n", config, sweep, e2);
    return sweep;
}

void init_gauge_transformation(mtrx_3x3 * restrict G){
    pos_vec position;
    #pragma omp parallel for num_threads(NUM_THREADS) private(position) schedule(dynamic)
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            position.t = t;
            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        set_identity_3x3(get_gaugetransf(G, position));
                    }
                }
            }
        }                
}