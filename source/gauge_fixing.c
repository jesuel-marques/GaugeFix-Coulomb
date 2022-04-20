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
#include "fourvector_field.h"  //	Calculation of A_mu(n) and related things
#include "math_ops.h"          //	Math operations

typedef enum {R, S, T} submatrix;

static void SU3_local_update_U(matrix_3x3_double *U, const pos_vec position, const matrix_3x3_double *g) {
    //	Updates U only at a given position

    matrix_3x3_double g_dagger;

    for (lorentz_index mu = 0; mu < DIM; mu++) {
        //	U'_mu(x)=g(x).U_mu(x).1 for red-black updates

        accumulate_left_product_3x3(g, get_link(U, position, mu));

        //	U'_mu(x-mu)=1.U_mu(x-mu).g_dagger(x) for red-black updates

        SU3_hermitean_conjugate(g, &g_dagger);

        accumulate_right_product_3x3(get_link(U, hop_position_negative(position, mu), mu), &g_dagger);
    }
}

static void SU3_calculate_w(matrix_3x3_double *U, const pos_vec position, matrix_3x3_double *w) {
    //	Calculates 	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision, following the notation in hep-lat/9306018
    //	returns result in w.

    set_to_null_3x3(w);  //	Initializing w(n)=0

    matrix_3x3_double u_dagger_rear;

    // w(n)	calculation

    for (lorentz_index mu = 0; mu < DIM - 1 ; mu++) {
        //	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision

        accumulate_3x3(get_link(U, position, mu), w);

        SU3_hermitean_conjugate(get_link(U, hop_position_negative(position, mu), mu), &u_dagger_rear);

        accumulate_3x3(&u_dagger_rear, w);
    }
}

unsigned calculate_sweeps_to_next_measurement(const double e2, const unsigned sweeps_to_measurement_e2) {
    return sweeps_to_measurement_e2 + (unsigned)(initial_sweeps_to_measurement_e2 * (1.0 - log10(e2) / log10(tolerance))) + 10;
}
    
static double SU3_calculate_e2(matrix_3x3_double *U) {
    //	Calculates e2 (defined in hep-lat/0301019v2),
    //	used to find out distance to the gauge-fixed situation.
 
    double e2 = 0.0;

    #pragma omp parallel for reduction (+:e2) num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            matrix_3x3_double div_A;
            matrix_SU3_alg div_A_components;
            pos_vec position;

            position.t = t;
            double e2_slice = 0.0;

            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        
                        SU3_divergence_A(U, position, &div_A);
                        decompose_algebra_SU3(&div_A, &div_A_components);

                        for (SU3_color_alg_index a = 1; a <= pow2(Nc)-1; a++) {
                            //	Normalized sum of the squares of the color components of the divergence of A.

                            e2_slice += pow2(div_A_components.m[a]);
                        }

                    }
                }
            }
            e2 += e2_slice;
        }

    e2 /= (VOLUME);

    return e2;
}

static void SU3_update_sub_LosAlamos(const matrix_3x3_double *matrix_SU3, submatrix sub, matrix_3x3_double *update_SU3) {
    SU3_color_index a, b;

    set_to_null_3x3(update_SU3);

    update_SU3 -> m[elm(2 - sub, 2 - sub)] = 1.0;
    
    a = sub == T ? 1 : 0;
    b = sub == R ? 1 : 2;

    matrix_2x2_ck matrix_SU2;

    matrix_SU2.m[0] =  (creal(matrix_SU3 -> m[elm(a, a)]) 
                      + creal(matrix_SU3 -> m[elm(b, b)]));
    matrix_SU2.m[1] = -(cimag(matrix_SU3 -> m[elm(a, b)]) 
                      + cimag(matrix_SU3 -> m[elm(b, a)]));
    matrix_SU2.m[2] = -(creal(matrix_SU3 -> m[elm(a, b)]) 
                      - creal(matrix_SU3 -> m[elm(b, a)]));
    matrix_SU2.m[3] = -(cimag(matrix_SU3 -> m[elm(a, a)]) 
                      - cimag(matrix_SU3 -> m[elm(b, b)]));

    SU2_projection(&matrix_SU2);

    update_SU3 -> m[elm(a, a)] =      matrix_SU2.m[0] 
                                + I * matrix_SU2.m[3];
    update_SU3 -> m[elm(a, b)] =      matrix_SU2.m[2]
                                + I * matrix_SU2.m[1];
    update_SU3 -> m[elm(b, a)] =     -matrix_SU2.m[2] 
                                + I * matrix_SU2.m[1];
    update_SU3 -> m[elm(b, b)] =      matrix_SU2.m[0] 
                                - I * matrix_SU2.m[3];
}

static void SU3_LosAlamos_common_block(const matrix_3x3_double *w, matrix_3x3_double *total_update) {
    //	Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
    //	algorithm for SU(3), with a division of the update matrix in submatrices
    //	following the Cabbibo-Marinari trick. Actual update is obtained after a number
    //	of "hits" to be performed one after another.

    matrix_3x3_double update;

    matrix_3x3_double updated_w;

    copy_3x3(w, &updated_w);
    set_to_identity_3x3(total_update);

    for (unsigned short hits = 1; hits <= maxhits; hits++) {
        //	Each hit contains the Cabbibo-Marinari subdivision
        for (submatrix sub = R; sub <= T; sub++) {
            //	Submatrices are indicated by numbers from 0 to 2

            SU3_update_sub_LosAlamos(&updated_w, sub, &update);
            accumulate_left_product_3x3(&update, &updated_w);

            accumulate_left_product_3x3(&update, total_update);
            //	Updates matrix to total_update. It is the
            //	accumulated updates from the hits.
        }
    }
}

static void SU3_gaugefixing_overrelaxation(matrix_3x3_double *U, const pos_vec position) {
    //	Generalization of the algorithm described in hep-lat/0301019v2, using the
    //	Cabbibo-Marinari submatrices trick.
    //	It updates the g at the given position.

    matrix_3x3_double w;

    SU3_calculate_w(U, position, &w);  //	Calculating w(n)=h(n) for red black subdivision

    matrix_3x3_double update_LA;

    SU3_LosAlamos_common_block(&w, &update_LA);

    /*	The above function determines update_LA which would be the naïve
   	update to bring the local function to its mininum. However,
 	using overrelaxation, which
   	means using update_LA^omega instead of  update_LA, where 1<omega<2 
    makes the converge faster. update_LA^omega is calculated using 
    the first two terms of the binomial expansion: 
    update_LA^omega=I+omega(update_LA-I)+...=I(1-omega)+omega*update_LA+...*/

    matrix_3x3_double update_OR;

    // update_OR = update_LA^omega = Proj_SU3((I(1-omega)+omega*update_LA)
    set_to_identity_3x3(&update_OR);
    substitution_multiplication_by_scalar_3x3(1.0 - omega_OR, &update_OR);
    substitution_multiplication_by_scalar_3x3(omega_OR, &update_LA);
    accumulate_3x3(&update_LA, &update_OR);

    projection_SU3(&update_OR);

    SU3_local_update_U(U, position, &update_OR);
}

unsigned SU3_gauge_fix(matrix_3x3_double *U, const unsigned short config) {
    //	Fix the gauge and follows the process by calculating e2;

    pos_vec position;

    double e2;

    unsigned sweep = 0;
    unsigned sweeps_to_measurement_e2 = initial_sweeps_to_measurement_e2;
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
                            !((position_is_even(position) + sweep) % 2) ?
                                //	Implementation of the checkerboard subdivision of the lattice
                                
                                SU3_gaugefixing_overrelaxation(U, position)
                                //  The actual gauge-fixing algorithm
                                                                    
                                : 0;
                            
                        }
                    }
                }
            }

        sweep++;

        !(sweep % (initial_sweeps_to_measurement_e2 / 5) ) ?
            printf("Sweeps in config %5d: %8d.\n", config, sweep) : 0; 

        if(sweep == sweeps_to_measurement_e2){
            e2 = SU3_calculate_e2(U);
	        printf("Sweeps in config %5d: %8d. e2: %3.2E \n", config, sweep, e2);            //	Gauge-fixing index, indicates how far we are to the Landau-gauge.
            //  It will be less than the tolerance,
            //	when the gauge-fixing is considered to be attained.
            //	Following the notation of hep-lat/0301019v2
            
            if (e2 <= tolerance) {
                break;

            } else {
                sweeps_to_measurement_e2 = calculate_sweeps_to_next_measurement(e2, sweeps_to_measurement_e2);
                //	No need to calculate e2 all the time
                //	because it will take some hundreds of sweeps
                //	to fix the gauge.
            }
        }

        !(sweep % sweeps_to_reunitarization) ? SU3_reunitarize(U) : 0;
                                               
    }
    SU3_reunitarize(U);

    printf("Sweeps needed to gauge-fix config %d: %d. e2: %3.2E \n", config, sweep, e2);
    return sweep;
}
