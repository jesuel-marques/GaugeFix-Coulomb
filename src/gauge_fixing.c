#include <stdio.h>
#include <tgmath.h>

#include <omp.h>

#include <fields.h>
#include <flags.h>
#include <four_potential.h>
#include <gauge_fixing.h>
#include <lattice.h>
#include <misc.h>
#include <SU3_ops.h>

#include <types.h>

//  Initial amount of sweeps to actually measure e2
#define ESTIMATE_SWEEPS_TO_CHECK_e2 1000
//  Maximum number of sweeps that program will perform to fix
#define MAX_SWEEPS_TO_FIX 100 * ESTIMATE_SWEEPS_TO_CHECK_e2
//  Amount of sweeps to reunitarize fields
#define SWEEPS_TO_REUNITARIZATION 250

//  Iterations hits for the maximization of Tr[w(n)]
#define MAX_HITS 2  

//	Implementation of the checkerboard subdivision of the lattice.
#define BLACKRED_SWEEP(position, sweep) ((POSITION_IS_EVEN(position)) ^ ((sweep) & 1))


bool validGaugeFixingParametersQ(GaugeFixingParameters gfix_param) {
    if(gfix_param.omega_OR  < 1 || gfix_param.omega_OR >= 2 ||
       gfix_param.tolerance < 0 || gfix_param.error    == 1   ) {

		return false;

	}
    return true;
}


GaugeFixingParameters initGaugeFixingParameters(const double tolerance, 
                                                  const double omega_OR) {
    GaugeFixingParameters gfix_param = {.tolerance = 0.0, 
                                        .omega_OR  = 0.0, 
                                        .error     = 0};

	gfix_param.tolerance 	= tolerance;
	gfix_param.omega_OR 	= omega_OR;

    if(!validGaugeFixingParametersQ(gfix_param)) {

        gfix_param.error = 1;

    }
	return gfix_param;
}


static void updateLocalUG(Mtrx3x3 * restrict U, 
                             Mtrx3x3 * restrict G, 
                             const PosVec position, 
                             const Mtrx3x3 * restrict update) {
    /* Updates U and G only at a given position */
    accumLeftProd3x3(update, getGaugetransf(G, position));

    Mtrx3x3 update_dagger;
    hermConj3x3(update, &update_dagger);

    LorentzIdx mu;
    LOOP_LORENTZ(mu) {

        //	U'_mu(x)=g(x).U_mu(x).1 for red-black updates

        accumLeftProd3x3(update, getLink(U, position, mu));

        //	U'_mu(x-mu)=1.U_mu(x-mu).g_dagger(x) for red-black updates

        accumRightProd3x3(getLink(U, hopPosMinus(position, mu), mu), 
                             &update_dagger);

    }
}


void updateGlobalU(Mtrx3x3 * restrict U, 
                      Mtrx3x3 * restrict G) {
    /* Updates U on the whole lattice*/
    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t) {
        Mtrx3x3 * g;
        Mtrx3x3 * u;
        Mtrx3x3 u_updated;

        PosVec position;

        position.t = t;

        LOOP_SPATIAL(position) {

            g = getGaugetransf(G, position);
            LorentzIdx mu;
            LOOP_LORENTZ(mu) {

                //	U'_mu(x)=g(x).U_mu(x).gdagger(x+mu)
                u = getLink(U, position, mu);
                
                prod_vuwdagger3x3(g, u, getGaugetransf(G, hopPosPlus(position, mu)),
                                  &u_updated);
                
                copy3x3(&u_updated, u);

            }
        }
    }
}


static inline void accumulateFrontHearLink3x3(Mtrx3x3 * restrict U, 
                                              const PosVec position, 
                                              LorentzIdx mu, 
                                              Mtrx3x3 * restrict w) {
    /* Accumulates the value of U_mu(n) and U_-mu(n) into w */
    Mtrx3x3 * u_front = getLink(U,             position,      mu);
    Mtrx3x3 * u_rear  = getLink(U, hopPosMinus(position, mu), mu);

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {
        
        w -> m[ELEM_3X3(a, b)] +=       u_front -> m[ELEM_3X3(a, b)]
                                 + conj(u_rear  -> m[ELEM_3X3(b, a)]);        
    
    }
}


inline static void calculate_w(Mtrx3x3 * restrict U, 
                               const PosVec position,
                               Mtrx3x3 * restrict w) {
    /* 	Calculates 	w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 
        for red black subdivision, following the notation in hep-lat/9306018
    	returns result in w.*/
    setNull3x3(w);  //	Initializing w(n)=0

    // w(n)	calculation

    LorentzIdx mu;
    LOOP_LORENTZ_SPATIAL(mu) {

        /* w(n) = sum_mu U_mu(n).1+U_dagger_mu(n-mu_hat).1 for red black subdivision */
        accumulateFrontHearLink3x3(U, position, mu, w);

    }
}


static int sweepsToNextMeasurement(double e2, double tolerance) { 
    return 10 + (unsigned)(ESTIMATE_SWEEPS_TO_CHECK_e2 * 
                        (1.0 - log10(e2) / log10(tolerance))); 
}


double calculateF(Mtrx3x3 * restrict U) {
    Mtrx3x3 U_acc;
    setNull3x3(&U_acc);

    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t) {

        PosVec position;
        position.t = t;

        LOOP_SPATIAL(position) {
            LorentzIdx mu;
            LOOP_LORENTZ_SPATIAL(mu) {

                accumulate3x3(getLink(U, position, mu), &U_acc);

            }
        }
    }

    return creal(trace3x3(&U_acc)) / (lattice_param.volume * Nc * (DIM - 1));
}


double calculateTheta(Mtrx3x3 * restrict U) { 
    Mtrx3x3 div_A, div_A_dagger, prod;

    double theta = 0.0;

    PosIndex t;
    LOOP_TEMPORAL_PARALLEL(t) {

        PosVec position;
        position.t = t;
        LOOP_SPATIAL(position) {     

            divergenceA(U, position, &div_A);
            hermConj3x3(&div_A, &div_A_dagger);
            prod3x3(&div_A, &div_A_dagger, &prod);
            theta += trace3x3(&prod);

        }
    }

    return creal(theta) / (double) (Nc * lattice_param.volume);
}

    
double calculate_e2(Mtrx3x3 * restrict U) {
    /* 	Calculates e2 (defined in hep-lat/0301019v2),
    	used to find out distance to the gauge-fixed situation. */
    double e2 = 0.0;
    PosIndex t;
    #pragma omp parallel for reduction (+:e2) num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        LOOP_TEMPORAL(t) {
            Mtrx3x3 div_A;
            MtrxSU3Alg div_A_components;

            double e2_slice = 0.0;

            PosVec position;

            position.t = t;

            LOOP_SPATIAL(position) {

                divergenceA(U, position, &div_A);
                decomposeAlgebraSU3(&div_A, &div_A_components);
                for (SU3AlgIdx a = 1; a <= POW2(Nc)-1; a++) {

                    /* 	Normalized sum of the squares of the color components
                        of the divergence of A. */
                    e2_slice += POW2(div_A_components.m[a]);

                }

            }
              
            e2 += e2_slice;
        }

    e2 /= (lattice_param.volume);

    return e2;
}


inline void updateSubLosAlamos(Mtrx3x3 * restrict w, 
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

    if(!projectSU2(&mtrx_SU2)) {

        /* If projectSU2 is succesful, it will return 0
           Then w can be updated */
        accumProdSU2_3x3(&mtrx_SU2, w, a, b);

    }
    else{
        /*  If projectSU2 is unsuccesful, this means
            that mtrx_SU2 was 0 and w remains what it was.
            The program will then carry on to the next 
            Submtrx update. */
        fprintf(stderr, "A matrix could not be projected to SU(2)");
        return;
    }
}


inline void updateLosAlamos(Mtrx3x3 * restrict w, 
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
    

       
    if(inverse3x3(w, &w_inv_old)) {

        /*  Local maximization is attained iteratively in SU(3),
            thus we need to make many hits ... */
        for (unsigned short hits = 1; hits <= MAX_HITS; hits++) {

            /* ... and each hit contains the Cabbibo-Marinari subdivision */
            for (Submtrx sub = R; sub <= T; sub++) {
                /* Submatrices are indicated by numbers from 0 to 2
                   with codenames R, S and T */

                updateSubLosAlamos(w, sub);
                
            }
        }
    }
    else{

        /* if w has no inverse, update will be given by the identity */
        setIdentity3x3(total_update);

    }
    
    prod3x3(w, &w_inv_old, total_update);
    /*	Updates matrix to total_update. It is the
    	accumulated updates from the hits. */
}


inline static void gaugefixOverrelaxationLocal(      Mtrx3x3 * restrict U,  
                                                     Mtrx3x3 * restrict G, 
                                                     const PosVec position, 
                                                     const double omega_OR) {
    /*	Generalization of the algorithm described in hep-lat/0301019v2, using the
    	Cabbibo-Marinari submatrices trick.
    	It updates the g at the given position.*/
    Mtrx3x3 w;

    //	Calculating w(n)=h(n) for red black subdivision
    calculate_w(U, position, &w);  

    Mtrx3x3 update_LA;

    updateLosAlamos(&w, &update_LA);

    /*	The above function determines update_LA which would be the naïve
   	update to bring the local function to its mininum. However,
 	using overrelaxation, which
   	means using update_LA^omega instead of  update_LA, where 1<omega<2 
    makes the converge faster. update_LA^omega is calculated using 
    the first two terms of the binomial expansion: 
    update_LA^omega=I+omega(update_LA-I)+...=I(1-omega)+omega*update_LA+...*/

    Mtrx3x3 update_OR;

    /* update_OR = update_LA^omega = Proj_SU3((I(1-omega)+omega*update_LA) */
    
    if(power3x3Binomial(&update_LA, omega_OR, &update_OR))
        setIdentity3x3(&update_OR);   
        /*  if matrix could not be projected to SU3 inside
        power_3x3 binomial, then use identity as update */

    updateLocalUG(U, G, position, &update_OR);
}


int gaugefixOverrelaxation(Mtrx3x3 * restrict U, 
                                Mtrx3x3 * restrict G, 
                                const GaugeFixingParameters gfix_param) {

    if(!validGaugeFixingParametersQ(gfix_param)) {

        return -2;

    }

    double last_e2 = 10.0;
    double new_e2 = calculate_e2(U);
    
    //	Counter to the number of sweeps to fix config to Coulomb gauge
    int sweep = 0;
    int sweeps_to_measurement_e2 = sweepsToNextMeasurement(new_e2, 
                                                           gfix_param.tolerance);

    while (new_e2 > gfix_param.tolerance) {
        OMP_PARALLEL_FOR
            // Paralelizing by slicing the time extent
            for (PosIndex t = 0; t < lattice_param.n_T; t++) {
                PosVec position;
                position.t = t;
                LOOP_SPATIAL(position) {

                    BLACKRED_SWEEP(position, sweep) ?
                        //	Implementation of the checkerboard subdivision 
                        //  of the lattice.
                        
                        gaugefixOverrelaxationLocal(U, G, position, 
                                                          gfix_param.omega_OR): 0;
                        //  The actual gauge-fixing algorithm

                }
            }

        sweep++;

        if(sweep >= MAX_SWEEPS_TO_FIX ) {

                return -1;

        }

        if(sweep == sweeps_to_measurement_e2) {

            last_e2 = new_e2;
            new_e2  = calculate_e2(U);

            if (new_e2 <= gfix_param.tolerance) {

                break;

            }
            
            else if(last_e2 == new_e2 ) {

                return -1;

            }
            else {

                sweeps_to_measurement_e2 += 
                    sweepsToNextMeasurement(new_e2, gfix_param.tolerance);
                //	No need to calculate e2 all the time
                //	because it will take some hundreds/thousands of sweeps
                //	to fix the gauge.
                
            }
            printf("Sweeps: %8d.\t e2: %3.2E \n",
                   sweep,
                   new_e2);            
            //	Gauge-fixing index, indicates how far we are to the Landau-gauge.
            //  It will be less than the tolerance,
            //	when the gauge-fixing is considered to be attained.
            //	Following the notation of hep-lat/0301019v2

        }

        if(!(sweep % SWEEPS_TO_REUNITARIZATION)) {

            reunitarizeField(U, lattice_param.amount_of_links);
            reunitarizeField(G, lattice_param.amount_of_points);

        }
    }

    reunitarizeField(U, lattice_param.amount_of_links);
    reunitarizeField(G, lattice_param.amount_of_points);

    return sweep;
}