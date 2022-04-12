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

static void SU3_local_update_U(double complex *U, const pos_vec position, const double complex *g) {
    //	Updates U only at a given position

    double complex g_dagger[3 * 3];

    for (unsigned short mu = 0; mu < d; mu++) {
        //	U'_mu(x)=g(x).U_mu(x).1 for red-black updates

        SU3_accumulate_left_product(g, get_link(U, position, mu));

        //	U'_mu(x-mu)=1.U_mu(x-mu).g_dagger(x) for red-black updates

        SU3_hermitean_conjugate(g, g_dagger);

        SU3_accumulate_right_product(get_link(U, hop_position_negative(position, mu), mu), g_dagger);
    }
}

static void SU3_calculate_w(double complex *U, const pos_vec position, double complex *w) {
    //	Calculates 	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision, following the notation in hep-lat/9306018
    //	returns result in w.

    SU3_set_to_null(w);  //	Initializing w(n)=0

    double complex u_dagger_rear[3 * 3];

    // w(n)	calculation

    for (unsigned short mu = 1; mu < d; mu++) {
        //	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision

        SU3_accumulate(get_link(U, position, mu), w);

        SU3_hermitean_conjugate(get_link(U, hop_position_negative(position, mu), mu), u_dagger_rear);

        SU3_accumulate(u_dagger_rear, w);
    }
}

unsigned calculate_sweeps_to_next_measurement(const double e2, const unsigned sweeps_to_measurement_e2) {
    return sweeps_to_measurement_e2 + (unsigned)(initial_sweeps_to_measurement_e2 * (1.0 - log10(e2) / log10(tolerance))) + 10;
}
    
static double SU3_calculate_e2(double complex *U) {
    //	Calculates e2 (defined in hep-lat/0301019v2),
    //	used to find out distance to the gauge-fixed situation.
 
    double e2 = 0.0;

    #pragma omp parallel for reduction (+:e2) num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (unsigned short t = 0; t < Nt; t++) {
            double complex div_A[3 * 3];
            double div_A_components[9];
            pos_vec position;

            position.t = t;
            double e2_slice = 0.0;

            for (position.i = 0; position.i < Nxyz; position.i++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.k = 0; position.k < Nxyz; position.k++) {
                        
                        SU3_divergence_A(U, position, div_A);
                        SU3_decompose_algebra(div_A, div_A_components);

                        for (unsigned short a = 1; a <= 8; a++) {
                            //	Normalized sum of the squares of the color components of the divergence of A.

                            e2_slice += (double)pow2(div_A_components[a]);
                        }

                    }
                }
            }
            e2 += e2_slice;
        }

    e2 /= (Volume);

    return e2;
}

static void SU3_update_sub_LosAlamos(const double complex *matrix_SU3, unsigned short submatrix, complex double *update_SU3) {
    unsigned short a, b;

    SU3_set_to_null(update_SU3);

    update_SU3[(2 - submatrix) * 3 + (2 - submatrix)] = 1.0;
    a = submatrix == 2 ? 1 : 0;
    b = submatrix == 0 ? 1 : 2;

    double matrix_SU2[4];

    matrix_SU2[0] =  (creal(matrix_SU3[a * 3 + a]) + creal(matrix_SU3[b * 3 + b]));
    matrix_SU2[1] = -(cimag(matrix_SU3[a * 3 + b]) + cimag(matrix_SU3[b * 3 + a]));
    matrix_SU2[2] = -(creal(matrix_SU3[a * 3 + b]) - creal(matrix_SU3[b * 3 + a]));
    matrix_SU2[3] = -(cimag(matrix_SU3[a * 3 + a]) - cimag(matrix_SU3[b * 3 + b]));

    SU2_projection(matrix_SU2);

    update_SU3[a * 3 + a] =  matrix_SU2[0] + I * matrix_SU2[3];
    update_SU3[a * 3 + b] =  matrix_SU2[2] + I * matrix_SU2[1];
    update_SU3[b * 3 + a] = -matrix_SU2[2] + I * matrix_SU2[1];
    update_SU3[b * 3 + b] =  matrix_SU2[0] - I * matrix_SU2[3];
}

static void SU3_LosAlamos_common_block(const double complex *w, double complex *total_update) {
    //	Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
    //	algorithm for SU(3), with a division of the update matrix in submatrices
    //	following the Cabbibo-Marinari trick. Actual update is obtained after a number
    //	of "hits" to be performed one after another.

    double complex update[3 * 3];

    double complex updated_w[3 * 3];

    SU3_copy(w, updated_w);
    SU3_set_to_identity(total_update);

    for (unsigned short hits = 1; hits <= maxhits; hits++) {
        //	Each hit contains the Cabbibo-Marinari subdivision
        for (unsigned short submatrix = 0; submatrix <= 2; submatrix++) {
            //	Submatrices are indicated by numbers from 0 to 2

            SU3_update_sub_LosAlamos(updated_w, submatrix, update);
            SU3_accumulate_left_product(update, updated_w);

            SU3_accumulate_left_product(update, total_update);
            //	Updates matrix to total_update. It is the
            //	accumulated updates from the hits.
        }
    }
}

static void SU3_gaugefixing_overrelaxation(double complex *U, const pos_vec position) {
    //	Generalization of the algorithm described in hep-lat/0301019v2, using the
    //	Cabbibo-Marinari submatrices trick.
    //	It updates the g at the given position.

    double complex w[3 * 3];

    SU3_calculate_w(U, position, w);  //	Calculating w(n)=h(n) for red black subdivision

    double complex update_LA[3 * 3];

    SU3_LosAlamos_common_block(w, update_LA);

    /*	The above function determines update_LA which would be the naïve
   	update to bring the local function to its mininum. However,
 	using overrelaxation, which
   	means using update_LA^omega instead of  update_LA, where 1<omega<2 
    makes the converge faster. update_LA^omega is calculated using 
    the first two terms of the binomial expansion: 
    update_LA^omega=I+omega(update_LA-I)+...=I(1-omega)+omega*update_LA+...*/

    double complex update_OR[3 * 3];

    // update_OR = update_LA^omega = Proj_SU3((I(1-omega)+omega*update_LA)
    SU3_set_to_identity(update_OR);
    SU3_substitution_multiplication_by_scalar(1.0 - omega_OR, update_OR);
    SU3_substitution_multiplication_by_scalar(omega_OR, update_LA);
    SU3_accumulate(update_LA, update_OR);

    SU3_projection(update_OR);

    SU3_local_update_U(U, position, update_OR);
}

unsigned SU3_gauge_fix(double complex *U, const unsigned short config) {
    //	Fix the gauge and follows the process by calculating e2;

    pos_vec position;

    double e2;

    unsigned sweep = 0;
    unsigned sweeps_to_measurement_e2 = initial_sweeps_to_measurement_e2;
    //	Counter to the number of sweeps to fix config to Landau gauge
	e2 = SU3_calculate_e2(U);
	printf("Sweeps in config %d: %d. e2: %3.2E \n", config, sweep, e2);
    
    while (1) {
        #pragma omp parallel for num_threads(NUM_THREADS) private(position) schedule(dynamic)
            // Paralelizing by slicing the time extent
            for (unsigned short t = 0; t < Nt; t++) {
                position.t = t;
                for (position.i = 0; position.i < Nxyz; position.i++) {
                    for (position.j = 0; position.j < Nxyz; position.j++) {
                        for (position.k = 0; position.k < Nxyz; position.k++) {
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

        !(sweep % 50) ?
            printf("sweep: %d\n", sweep) : 0; 

        if(sweep == sweeps_to_measurement_e2){
            e2 = SU3_calculate_e2(U);
            printf("Sweeps in config %d: %d. e2: %3.2E \n", config, sweep, e2);
            //	Gauge-fixing index, indicates how far we are to the Landau-gauge.
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
