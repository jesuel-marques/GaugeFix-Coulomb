#include <stdio.h>  //	Standard header C files
#include <stdlib.h>
#include <tgmath.h>

#include <omp.h>

#include <SU3_parameters.h>              //	Simulation parameters
#include <SU3_gaugefixing_parameters.h>  //	Gauge-fixing specific parameters
#include <fields.h>                        //	Initialization functions and
                                            //	calculations of positions and
                                            //	links on the lattice.

#include <math_ops.h>          //	Math operations
#include <SU2_ops.h>           //	SU(2) operations
#include <SU3_ops.h>           //	SU(3) operations
#include <four_potential.h>  //	Calculation of A_mu(n) and related things
#include <gauge_fixing.h>
#include <misc.h>
#include <lattice.h>


#define SWEEPS_TO_NEXT_MEASUREMENT(e2)  10 + (unsigned)(INITIAL_SWEEPS_TO_MEASUREMENT_e2 \
                                              * (1.0 - log10((e2)) / log10(TOLERANCE)))

#define CHECKERBOARD_SUBDIVISION(position, sweep) !((POSITION_IS_EVEN(position)) ^ ((sweep) & 1))
//	Implementation of the checkerboard subdivision of the lattice.

static void SU3_local_update_U_G(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G, const pos_vec position, 
                                                        const mtrx_3x3 * restrict update) {
    //	Updates U and G only at a given position

    accum_left_prod_3x3(update, get_gaugetransf(G, position));

    mtrx_3x3 update_dagger;
    herm_conj_3x3(update, &update_dagger);

    for (lorentz_idx mu = 0; mu < DIM; mu++) {
        //	U'_mu(x)=g(x).U_mu(x).1 for red-black updates

        accum_left_prod_3x3(update, get_link(U, position, mu));

        //	U'_mu(x-mu)=1.U_mu(x-mu).g_dagger(x) for red-black updates

        accum_right_prod_3x3(get_link(U, hop_pos_minus(position, mu), mu), 
                                                                        &update_dagger);
    }
}

void SU3_global_update_U(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G) {
    //	Updates U on the whole lattice

    OMP_PARALLEL_FOR
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            mtrx_3x3 * g;
            mtrx_3x3 * u;
            mtrx_3x3 u_updated;

            pos_vec position;

            position.t = t;

            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        g = get_gaugetransf(G, position);
                        for (lorentz_idx mu = 0; mu < DIM; mu++) {
                            //	U'_mu(x)=g(x).U_mu(x).gdagger(x+mu)
                            u = get_link(U, position, mu);
                            
                            prod_vuwdagger_3x3(g, u, get_gaugetransf(G, hop_pos_plus(position, mu)), &u_updated);
                            
                            copy_3x3(&u_updated, u);

                        }
                    }
                }
            }
        }
}


static inline void accumulate_front_hear_link_3x3(mtrx_3x3 * restrict U, const pos_vec position, lorentz_idx mu, 
                                                                                        mtrx_3x3 * restrict w) {
    // Accumulates the value of U_mu(n) and U_-mu(n) into w

    mtrx_3x3 * u_front = get_link(U,               position,      mu);
    mtrx_3x3 * u_rear  = get_link(U, hop_pos_minus(position, mu), mu);

    SU3_color_idx a, b;
    LOOP_COLOR_3x3(a, b){

            w -> m[ELM3x3(a, b)] +=     u_front -> m[ELM3x3(a, b)]
                                 + conj(u_rear  -> m[ELM3x3(b, a)]);
        
        
    }
}


inline static void SU3_calculate_w(mtrx_3x3 * restrict U, const pos_vec position,
                                                         mtrx_3x3 * restrict w) {
    //	Calculates 	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision, following the notation in hep-lat/9306018
    //	returns result in w.

    set_null_3x3(w);  //	Initializing w(n)=0

    // w(n)	calculation

    for (lorentz_idx mu = 0; mu < DIM - 1 ; mu++) {
        //	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision

        accumulate_front_hear_link_3x3(U, position, mu, w);

    }
}

double SU3_calculate_F(mtrx_3x3 * restrict U){

    pos_vec position;

    mtrx_3x3 U_acc;
    set_null_3x3(&U_acc);
    
    for (position.t = 0; position.t < N_T; position.t++)  {
        for (position.k = 0; position.k < N_SPC; position.k++) {
            for (position.j = 0; position.j < N_SPC; position.j++) {
                for (position.i = 0; position.i < N_SPC; position.i++) {
                    for (lorentz_idx mu = 0; mu < DIM - 1 ; mu++) {                        
                        accumulate_3x3(get_link(U, position, mu), &U_acc);                    
                    }
                }
            }
        }
    }

    return creal(trace_3x3(&U_acc)) / (VOLUME * Nc * (DIM - 1));
}

double SU3_calculate_theta(mtrx_3x3 * restrict U){
    
    pos_vec position;

    mtrx_3x3 div_A, div_A_dagger, prod;

    double theta = 0.0;

    for (position.t = 0; position.t < N_T; position.t++)  {
        for (position.k = 0; position.k < N_SPC; position.k++) {
            for (position.j = 0; position.j < N_SPC; position.j++) {
                for (position.i = 0; position.i < N_SPC; position.i++) {

                    SU3_divergence_A(U, position, &div_A);
                    herm_conj_3x3(&div_A, &div_A_dagger);
                    prod_3x3(&div_A, &div_A_dagger, &prod);
                    theta += trace_3x3(&prod);

                }
            }
        }
    }

    return creal(theta) / (Nc * VOLUME);
}
    
double SU3_calculate_e2(mtrx_3x3 * restrict U) {
    //	Calculates e2 (defined in hep-lat/0301019v2),
    //	used to find out distance to the gauge-fixed situation.
 
    double e2 = 0.0;

    #pragma omp parallel for reduction (+:e2) num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            mtrx_3x3 div_A;
            mtrx_SU3_alg div_A_components;

            double e2_slice = 0.0;

            pos_vec position;

            position.t = t;

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

inline void SU3_update_sub_LosAlamos(mtrx_3x3 * restrict w, submatrix sub) {
    
    SU3_color_idx a, b;

    a = sub == T ? 1 : 0;
    b = sub == R ? 1 : 2;

    //  a and b will be the line and colums index
    //  for each Cabbibo-Marinari matrix, 

    mtrx_2x2_ck mtrx_SU2;

    mtrx_SU2.m[0] =  creal( w -> m[ELM3x3(a, a)]
                          + w -> m[ELM3x3(b, b)]);
    mtrx_SU2.m[1] = -cimag( w -> m[ELM3x3(a, b)] 
                          + w -> m[ELM3x3(b, a)]);
    mtrx_SU2.m[2] = -creal( w -> m[ELM3x3(a, b)]
                          - w -> m[ELM3x3(b, a)]);
    mtrx_SU2.m[3] = -cimag( w -> m[ELM3x3(a, a)] 
                          - w -> m[ELM3x3(b, b)]);
  

    //  The SU(2) matrix corresponding to the Cabbibo-Marinari 
    //  submatrix is built from w according to the formulae above.
    //  This is what maximizes the functional locally for each
    //  submatrix. This can be proven using maximization with 
    //  constrains, where the constraint is that mtrx_SU2 has 
    //  to be an SU(2) matrix.

    if(!SU2_projection(&mtrx_SU2)){
        //  If SU2_projection is succesful, it will return 0
        //  Then w can be updated
        accum_prod_SU2_3x3(&mtrx_SU2, w, a, b);
    }
    else{
        //  If SU2_projection is unsuccesful, this means
        //  that mtrx_SU2 was 0 and w remains what it was.
        //  The program will then carry on to the next 
        //  submatrix update.
        return;
    }
}

inline void SU3_LosAlamos_common_block(mtrx_3x3 * restrict w, 
                                              mtrx_3x3 * restrict total_update) {
    //	Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
    //	algorithm for SU(3), with a division of the update matrix in submatrices
    //	following the Cabbibo-Marinari trick. Actual update is obtained after a number
    //	of "hits" to be performed one after another.

    mtrx_3x3 w_inv_old; 
    //  Calculates the inverse of w in the beginning.
    //  The program will update w successively and to
    //  extract what was the combined update, we can 
    //  multiply from the right by the old inverse.
    

       
    if(inverse_3x3(w, &w_inv_old)){

        

        //  Local maximization is attained iteratively in SU(3),
        //  thus we need to make many hits ...
        for (unsigned short hits = 1; hits <= MAX_HITS; hits++) {

            //	... and each hit contains the Cabbibo-Marinari subdivision
            for (submatrix sub = R; sub <= T; sub++) {
                //	Submatrices are indicated by numbers from 0 to 2
                //  with codenames R, S and T

                SU3_update_sub_LosAlamos(w, sub);
                
            }
        }
    }
    else{
        //  if w has no inverse, update will be given by the identity
        set_identity_3x3(total_update);
    }
    
    prod_3x3(w, &w_inv_old, total_update);
    //	Updates matrix to total_update. It is the
    //	accumulated updates from the hits.
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
    
    if(power_3x3_binomial(&update_LA, OMEGA_OR, &update_OR))
        set_identity_3x3(&update_OR);   
        /*  if matrix could not be projected to SU3 inside
        power_3x3 binomial, then use identity as update */

    SU3_local_update_U_G(U, G, position, &update_OR);
}

int SU3_gauge_fix(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G, const unsigned short config_nr) {
    //	Fix the gauge and follows the process by calculating e2;


    double last_e2 = 10, new_e2 = 10;

    int sweep = 0;
    int sweeps_to_measurement_e2 = INITIAL_SWEEPS_TO_MEASUREMENT_e2;
    //	Counter to the number of sweeps to fix config to Landau gauge

	new_e2 = SU3_calculate_e2(U);
	printf("Sweeps in config %5d: %8d. e2: %3.2E \n", config_nr, sweep, new_e2);

    while (1) {

        OMP_PARALLEL_FOR
            // Paralelizing by slicing the time extent
            for (pos_index t = 0; t < N_T; t++) {
                pos_vec position;
                position.t = t;
                for (position.k = 0; position.k < N_SPC; position.k++) {
                    for (position.j = 0; position.j < N_SPC; position.j++) {
                        for (position.i = 0; position.i < N_SPC; position.i++) {
                            CHECKERBOARD_SUBDIVISION(position, sweep) ?
                                //	Implementation of the checkerboard subdivision 
                                //  of the lattice.
                                
                                SU3_gaugefixing_overrelaxation(U, G, position)
                                //  The actual gauge-fixing algorithm
                                                                    
                                : 0;
                            
                        }
                    }
                }
            }

        // printf("sweep: %d\n", sweep);

        sweep++;

        // if(!system("test -f stop_run")){

        //     printf("Request to stop run for config %d in sweep %d.\n", config_nr, sweep);
        //     return -1;
        
        // }


        if(sweep == sweeps_to_measurement_e2) {
            new_e2 = SU3_calculate_e2(U);

            if (new_e2 <= TOLERANCE) {

                break;

            }
            
            else if(sweep >= 100 * INITIAL_SWEEPS_TO_MEASUREMENT_e2 || last_e2 == new_e2 ){

                fprintf(stderr, "Config could not be gauge-fixed.\n SOR algorithm seems not to work or be too slow for config %d.\n", config_nr);
                return -1;

            }
            else {
                
                sweeps_to_measurement_e2 += 
                    SWEEPS_TO_NEXT_MEASUREMENT(new_e2);
                //	No need to calculate e2 all the time
                //	because it will take some hundreds of sweeps
                //	to fix the gauge.

            }
            printf("Sweeps in config %5d: %8d. e2: %3.2E \n", config_nr, sweep, new_e2);            
            //	Gauge-fixing index, indicates how far we are to the Landau-gauge.
            //  It will be less than the TOLERANCE,
            //	when the gauge-fixing is considered to be attained.
            //	Following the notation of hep-lat/0301019v2
            
        }
        else if(!(sweep % (INITIAL_SWEEPS_TO_MEASUREMENT_e2 / 2) )) {

            printf("Sweeps in config %5d: %8d.\n", config_nr, sweep);

        }

        (sweep % SWEEPS_TO_REUNITARIZATION) ? 0 : SU3_reunitarize_U_G(U, G);

        last_e2 = new_e2;
    }
    SU3_reunitarize_U_G(U, G);

    printf("Sweeps needed to gauge-fix config %d: %d. e2: %3.2E \n", config_nr, sweep, new_e2);
    return sweep;
}

void init_gauge_transformation(mtrx_3x3 * restrict G){
    
    OMP_PARALLEL_FOR
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {
            pos_vec position;

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