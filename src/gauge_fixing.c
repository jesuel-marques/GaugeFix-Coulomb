#include <stdio.h>
#include <tgmath.h>

#include <omp.h>

#include <fields.h>
#include <flags.h>
#include <four_potential.h>
#include <gauge_fixing.h>
#include <lattice.h>
#include <math_ops.h>         
#include <misc.h>
#include <SU3_ops.h>
#include <SU3_parameters.h>
#include <types.h>

//  Initial amount of sweeps to actually measure e2
#define INITIAL_SWEEPS_TO_CHECK_e2 1000
//  Amount of sweeps to reunitarize fields
#define SWEEPS_TO_REUNITARIZATION 250

//  Iterations hits for the maximization of Tr[w(n)]
#define MAX_HITS 2  

//	Implementation of the checkerboard subdivision of the lattice.
#define BLACKRED_SWEEP(position, sweep) ((POSITION_IS_EVEN(position)) ^ ((sweep) & 1))


extern int volume ;

extern short n_T;
extern short n_SPC;

extern int amount_of_links;
extern int amount_of_points;

static void SU3_local_update_U_G(Mtrx3x3 * restrict U, 
                                 Mtrx3x3 * restrict G, 
                                 const PosVec position, 
                                 const Mtrx3x3 * restrict update) {
    /* Updates U and G only at a given position */
    accum_left_prod_3x3(update, get_gaugetransf(G, position));

    Mtrx3x3 update_dagger;
    herm_conj_3x3(update, &update_dagger);

    LorentzIdx mu;
    LOOP_LORENTZ(mu) {
        //	U'_mu(x)=g(x).U_mu(x).1 for red-black updates

        accum_left_prod_3x3(update, get_link(U, position, mu));

        //	U'_mu(x-mu)=1.U_mu(x-mu).g_dagger(x) for red-black updates

        accum_right_prod_3x3(get_link(U, hop_pos_minus(position, mu), mu), 
                             &update_dagger);
    }
}


void SU3_global_update_U(Mtrx3x3 * restrict U, 
                         Mtrx3x3 * restrict G) {
    /* Updates U on the whole lattice*/
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t){
        Mtrx3x3 * g;
        Mtrx3x3 * u;
        Mtrx3x3 u_updated;

        PosVec position;

        position.t = t;

        LOOP_SPATIAL(position){
            g = get_gaugetransf(G, position);
            LorentzIdx mu;
            LOOP_LORENTZ(mu) {
                //	U'_mu(x)=g(x).U_mu(x).gdagger(x+mu)
                u = get_link(U, position, mu);
                
                prod_vuwdagger_3x3(g, u, get_gaugetransf(G, hop_pos_plus(position, mu)),
                                   &u_updated);
                
                copy_3x3(&u_updated, u);

            }
        }
    }
}


static inline void accumulate_front_hear_link_3x3(      Mtrx3x3 * restrict U, 
                                                  const PosVec position, 
                                                        LorentzIdx mu, 
                                                        Mtrx3x3 * restrict w) {
    /* Accumulates the value of U_mu(n) and U_-mu(n) into w */
    Mtrx3x3 * u_front = get_link(U,               position,      mu);
    Mtrx3x3 * u_rear  = get_link(U, hop_pos_minus(position, mu), mu);

    MtrxIdx3 a, b;
    LOOP_3X3(a, b){

        w -> m[ELEM_3X3(a, b)] +=       u_front -> m[ELEM_3X3(a, b)]
                                 + conj(u_rear  -> m[ELEM_3X3(b, a)]);
        
        
    }
}


inline static void SU3_calculate_w(      Mtrx3x3 * restrict U, 
                                   const PosVec position,
                                         Mtrx3x3 * restrict w) {
    /* 	Calculates 	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 
        for red black subdivision, following the notation in hep-lat/9306018
    	returns result in w.*/
    set_null_3x3(w);  //	Initializing w(n)=0

    // w(n)	calculation

    LorentzIdx mu;
    LOOP_LORENTZ_SPATIAL(mu){
        /* w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision */

        accumulate_front_hear_link_3x3(U, position, mu, w);

    }
}


static int sweeps_to_next_measurement(double e2, double tolerance){ 
    return 10 + (unsigned)(INITIAL_SWEEPS_TO_CHECK_e2 * 
                           (1.0 - log10(e2) / log10(tolerance))); 
}


double SU3_calculate_F(Mtrx3x3 * restrict U){
    Mtrx3x3 U_acc;
    set_null_3x3(&U_acc);

    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t){

        PosVec position;
        position.t = t;

        LOOP_SPATIAL(position){
            LorentzIdx mu;
            LOOP_LORENTZ_SPATIAL(mu){            
                accumulate_3x3(get_link(U, position, mu), &U_acc);                    
            }
        }
    }

    return creal(trace_3x3(&U_acc)) / (volume * Nc * (DIM - 1));
}


double SU3_calculate_theta(Mtrx3x3 * restrict U){ 
    Mtrx3x3 div_A, div_A_dagger, prod;

    double theta = 0.0;

    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t){

        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position){     

            SU3_divergence_A(U, position, &div_A);
            herm_conj_3x3(&div_A, &div_A_dagger);
            prod_3x3(&div_A, &div_A_dagger, &prod);
            theta += trace_3x3(&prod);

        }
    }

    return creal(theta) / (double) (Nc * volume);
}

    
double SU3_calculate_e2(Mtrx3x3 * restrict U) {
    /* 	Calculates e2 (defined in hep-lat/0301019v2),
    	used to find out distance to the gauge-fixed situation. */
    double e2 = 0.0;
    #pragma omp parallel for reduction (+:e2) num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (PosIndex t = 0; t < n_T; t++) {
            Mtrx3x3 div_A;
            MtrxSU3Alg div_A_components;

            double e2_slice = 0.0;

            PosVec position;

            position.t = t;

            LOOP_SPATIAL(position){

                SU3_divergence_A(U, position, &div_A);
                decompose_algebra_SU3(&div_A, &div_A_components);
                for (SU3AlgIdx a = 1; a <= POW2(Nc)-1; a++) {
                    /* 	Normalized sum of the squares of the color components
                        of the divergence of A. */
                    e2_slice += POW2(div_A_components.m[a]);
                }
                
            }
              
            e2 += e2_slice;
        }

    e2 /= (volume);

    return e2;
}


inline void SU3_update_sub_LosAlamos(Mtrx3x3 * restrict w, 
                                     Submtrx sub) {
    MtrxIdx3 a, b;

    a = sub == T ? 1 : 0;
    b = sub == R ? 1 : 2;

    /* a and b will be the line and colums index
    for each Cabbibo-Marinari matrix,  */

    Mtrx2x2CK mtrx_SU2;

    mtrx_SU2.m[0] =  creal( w -> m[ELEM_3X3(a, a)]
                          + w -> m[ELEM_3X3(b, b)]);
    mtrx_SU2.m[1] = -cimag( w -> m[ELEM_3X3(a, b)] 
                          + w -> m[ELEM_3X3(b, a)]);
    mtrx_SU2.m[2] = -creal( w -> m[ELEM_3X3(a, b)]
                          - w -> m[ELEM_3X3(b, a)]);
    mtrx_SU2.m[3] = -cimag( w -> m[ELEM_3X3(a, a)] 
                          - w -> m[ELEM_3X3(b, b)]);
  

    /*  The SU(2) matrix corresponding to the Cabbibo-Marinari 
        Submtrx is built from w according to the formulae above.
        This is what maximizes the functional locally for each
        Submtrx. This can be proven using maximization with 
        constrains, where the constraint is that mtrx_SU2 has 
        to be an SU(2) matrix. */

    if(!SU2_projection(&mtrx_SU2)){
        /* If SU2_projection is succesful, it will return 0
           Then w can be updated */
        accum_prod_SU2_3x3(&mtrx_SU2, w, a, b);
    }
    else{
        /*  If SU2_projection is unsuccesful, this means
            that mtrx_SU2 was 0 and w remains what it was.
            The program will then carry on to the next 
            Submtrx update. */
        fprintf(stderr, "A matrix could not be projected to SU(2)");
        return;
    }
}


inline void SU3_LosAlamos(Mtrx3x3 * restrict w, 
                          Mtrx3x3 * restrict total_update) {
    /*  Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
        algorithm for SU(3), with a division of the update matrix in submatrices
        following the Cabbibo-Marinari trick. Actual update is obtained after a number
        of "hits" to be performed one after another. */

    Mtrx3x3 w_inv_old; 
     /* First calculates the inverse of w.
     The program will update w successively and to
     extract what was the combined update, we can 
     multiply from the right by the old inverse. */
    

       
    if(inverse_3x3(w, &w_inv_old)){

        /*  Local maximization is attained iteratively in SU(3),
            thus we need to make many hits ... */
        for (unsigned short hits = 1; hits <= MAX_HITS; hits++) {

            /* ... and each hit contains the Cabbibo-Marinari subdivision */
            for (Submtrx sub = R; sub <= T; sub++) {
                /* Submatrices are indicated by numbers from 0 to 2
                   with codenames R, S and T */

                SU3_update_sub_LosAlamos(w, sub);
                
            }
        }
    }
    else{
        /* if w has no inverse, update will be given by the identity */
        set_identity_3x3(total_update);
    }
    
    prod_3x3(w, &w_inv_old, total_update);
    /*	Updates matrix to total_update. It is the
    	accumulated updates from the hits. */
}


inline static void SU3_gaugefix_overrelaxation_local(      Mtrx3x3 * restrict U,  
                                                           Mtrx3x3 * restrict G, 
                                                     const PosVec position, 
                                                     const double omega_OR) {
    /*	Generalization of the algorithm described in hep-lat/0301019v2, using the
    	Cabbibo-Marinari submatrices trick.
    	It updates the g at the given position.*/
    Mtrx3x3 w;

    //	Calculating w(n)=h(n) for red black subdivision
    SU3_calculate_w(U, position, &w);  

    Mtrx3x3 update_LA;

    SU3_LosAlamos(&w, &update_LA);

    /*	The above function determines update_LA which would be the naïve
   	update to bring the local function to its mininum. However,
 	using overrelaxation, which
   	means using update_LA^omega instead of  update_LA, where 1<omega<2 
    makes the converge faster. update_LA^omega is calculated using 
    the first two terms of the binomial expansion: 
    update_LA^omega=I+omega(update_LA-I)+...=I(1-omega)+omega*update_LA+...*/

    Mtrx3x3 update_OR;

    /* update_OR = update_LA^omega = Proj_SU3((I(1-omega)+omega*update_LA) */
    
    if(power_3x3_binomial(&update_LA, omega_OR, &update_OR))
        set_identity_3x3(&update_OR);   
        /*  if matrix could not be projected to SU3 inside
        power_3x3 binomial, then use identity as update */

    SU3_local_update_U_G(U, G, position, &update_OR);
}


int SU3_gaugefix_overrelaxation(Mtrx3x3 * restrict U, 
                                Mtrx3x3 * restrict G, 
                                char * config_filename, 
                                double tolerance, 
                                double omega_OR) {
    //	Fix the gauge and follows the process by calculating e2;
    double last_e2 = 10.0;
    double new_e2 = SU3_calculate_e2(U);
    
    int sweep = 0;
    int sweeps_to_measurement_e2 = INITIAL_SWEEPS_TO_CHECK_e2;
    //	Counter to the number of sweeps to fix config to Landau gauge


    while (1) {

        OMP_PARALLEL_FOR
            // Paralelizing by slicing the time extent
            for (PosIndex t = 0; t < n_T; t++) {
                PosVec position;
                position.t = t;
                LOOP_SPATIAL(position){
                    BLACKRED_SWEEP(position, sweep) ?
                        //	Implementation of the checkerboard subdivision 
                        //  of the lattice.
                        
                        SU3_gaugefix_overrelaxation_local(U, G, position, omega_OR): 0;
                        //  The actual gauge-fixing algorithm

                }
            }


        sweep++;


        if(sweep == sweeps_to_measurement_e2) {
            last_e2 = new_e2;
            new_e2  = SU3_calculate_e2(U);

            if (new_e2 <= tolerance) {

                break;

            }
            
            else if(sweep >= 100 * INITIAL_SWEEPS_TO_CHECK_e2 || last_e2 == new_e2 ){

                return -1;

            }
            else {
                
                sweeps_to_measurement_e2 += 
                    sweeps_to_next_measurement(new_e2, tolerance);
                //	No need to calculate e2 all the time
                //	because it will take some hundreds/thousands of sweeps
                //	to fix the gauge.

            }
            printf("Sweeps in config from file %s: %8d. e2: %3.2E \n", config_filename, 
                                                                       sweep,
                                                                       new_e2);            
            //	Gauge-fixing index, indicates how far we are to the Landau-gauge.
            //  It will be less than the tolerance,
            //	when the gauge-fixing is considered to be attained.
            //	Following the notation of hep-lat/0301019v2
            
        }
        else if(!(sweep % (INITIAL_SWEEPS_TO_CHECK_e2 / 2) )) {

            // printf("Sweeps in config %5d: %8d.\n", config_nr, sweep);

        }

        if(!(sweep % SWEEPS_TO_REUNITARIZATION)){
            reunitarize_field(U, amount_of_links);
            reunitarize_field(G, amount_of_points);

        }

    }
    reunitarize_field(U, amount_of_links);
    reunitarize_field(G, amount_of_points);

    return sweep;
}