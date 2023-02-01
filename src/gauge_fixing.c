/*
    gauge_fixing contain routines to gauge fix a configuration.

    Copyright (C) 2023  Jesuel Marques

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Contact: jesuel.leal@usp.br
    
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include <omp.h>

#include <fields.h>
#include <../settings.h>
#include <four_potential.h>
#include <gauge_fixing.h>
#include <geometry.h>
#include <SU3_ops.h>
#include <types.h>

#define PARAMETER_KEY_VALUE_DELIMITER ":"

//	Implementation of the checkerboard subdivision of the lattice.
#define BLACKRED_SWEEP(position, sweep) ((POSITION_IS_EVEN(position)) ^ ((sweep) & 1))

/* Verifies if gfix_param contain valid parameters for the gauge-fixing. */
bool validORGaugeFixingParametersQ(ORGaugeFixingParameters gfix_param) {

    /* 
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
	 *
	 * Parameters:
	 * ==========
	 * ORGaugeFixingParameters	gfix_param:		overrelaxation gauge fixing parameters.
	 *
	 * Returns:
	 * =======
	 * true if parameters are valid, false if invalid.
	 */

    return (gfix_param.omega_OR  >= 1 && gfix_param.omega_OR < 2    && 
            gfix_param.generic_gf.tolerance > 0                     && 
            gfix_param.generic_gf.error != true                     &&
            gfix_param.generic_gf.sweeps_to_reunitarization > 0         );
}

ORGaugeFixingParameters initParametersORDefault(){
    ORGaugeFixingParameters gfix_param = 
                    {.omega_OR = 1.95,  
                     .hits = 2,                     
                     .generic_gf.error = false,
                     .generic_gf.gfix_proxy = calculate_e2,
                     .generic_gf.tolerance = 1E-16,
                     .generic_gf.max_sweeps_to_fix = 100000,
                     .generic_gf.estimate_sweeps_to_gf_progress = 1000,
                     .generic_gf.sweeps_to_reunitarization = 250};

    return gfix_param;
}

void printORGaugeFixingParameters(ORGaugeFixingParameters gfix_param){

    printf("\n");
    printf("max_hits: %u\n",
            gfix_param.hits);
    printf("omega_OR: %lf\n",
            gfix_param.omega_OR);
    printf("\n");
    printf("tolerance: %3.2E\n",
            gfix_param.generic_gf.tolerance);
    printf("max_sweeps_to_fix: %u\n", 
            gfix_param.generic_gf.max_sweeps_to_fix);
    printf("estimate_sweeps_to_gfixprogress: %u\n", 
            gfix_param.generic_gf.estimate_sweeps_to_gf_progress);
    printf("\n");
    printf("sweeps_to_reunitarization: %u\n",
            gfix_param.generic_gf.sweeps_to_reunitarization);
    printf("\n");
}

/* Removes whitespace characters from a string. */
void removeSpaces(char * restrict string)
{
    /*
	 * Calls:
	 * =====
	 * isspace.
	 *
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
	 * char * string:      string whose spaces will be removed           
     * 
	 * Returns:
	 * =======
	 * 
	 */

    int non_space_count = 0;

    for(unsigned count = 0; string[count]; count++) {
        if(!isspace(string[count])) {
            string[non_space_count] = string[count];
            non_space_count++;
        }
    }
    
    string[non_space_count] = '\0';
}


/* Initializes gauge-fixing parameters given a tolerance and a omega parameter
   for the overrelaxation algorithm, checking if the parameters passed are valid. */
ORGaugeFixingParameters initORGaugeFixingParameters(const char * parameter_filename) {

    /*
	 * Calls:
	 * =====
	 * validORGaugeFixingParametersQ.
	 *
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
	 * double omega_OR:                 omega parameter of overrelaxation algorithm,
     * unsigned hits:                   iterations hits for the maximization of Tr[w],
     * double (*gfix_prox)(Mtrx3x3 * ): a function pointer to the function which 
     *                                  calculates the residue, which will serve as a
     *                                  proxy of the gauge-fixing condition,
     * double tolerance:                tolerance for the gauge-fixing residue,
     * unsigned max_sweeps_to_fix:      
     * unsigned estimate_sweeps_to_gf_progress:
     * 
     * Initial estimate of sweeps it will take for gauge-fixing to progress substantially
     *  Amount of sweeps to reunitarize fields
     * 
	 * Returns:
	 * =======
	 * ORGaugeFixingParameters struct with parameters set. If parameters are invalid, 
     * the error flag in the struct will be set.
	 */

    ORGaugeFixingParameters gfix_param = initParametersORDefault();

    if(!strcmp(parameter_filename, " ")){
        printf("No parameter-file detected, using default values.\n");
        return gfix_param;
    }

    FILE * parameter_file = fopen(parameter_filename, "r");

    char delimiter[] = PARAMETER_KEY_VALUE_DELIMITER;
    char buff[100];

    if(parameter_file != NULL) {
        while(fgets(buff, 100, parameter_file)) {

            removeSpaces(buff);

            char * key;
            key = strtok(buff, delimiter);
            
            char * value;
            value = strtok(NULL, "\n");
            if(value){
                
                if(!strcmp(key, "tolerance")) {
                    gfix_param.generic_gf.tolerance = atof(value);            
                }else if(!strcmp(key, "hits")) {
                    gfix_param.hits = atof(value);
                }else if(!strcmp(key, "omega_OR")) {
                    gfix_param.omega_OR = atof(value);
                }else if(!strcmp(key, "max_sweeps_to_fix")) {
                    gfix_param.generic_gf.max_sweeps_to_fix = atof(value);
                }else if(!strcmp(key, "estimate_sweeps_to_gf_progress")) {
                    gfix_param.generic_gf.estimate_sweeps_to_gf_progress = atof(value);
                }else if(!strcmp(key, "sweeps_to_reunitarization")) {
                    gfix_param.generic_gf.sweeps_to_reunitarization = atof(value);
                }else if(!strcmp(key, "gfix_proxy")) {
                
                    if(!strcmp(value, "e2")){
                        gfix_param.generic_gf.gfix_proxy = calculate_e2;
                    }
                    else if(!strcmp(value, "theta")) {
                        gfix_param.generic_gf.gfix_proxy = calculateTheta;
                    }
                    else{
                        fprintf(stderr, "Unrecognized gfix_proxy: %s.\n", value);
                    }

                }else if(strcmp(key, "\n")) {
                    fprintf(stderr, "Unrecognized parameter: %s.\n", key);
                }

            }
        }
        fclose(parameter_file);
    }
    else{
        fprintf(stderr, "Gauge-fixing parameter file %s could not be opened. "
                        "Using default parameters instead.\n", parameter_filename);
        printORGaugeFixingParameters(gfix_param);
    }

    if(!validORGaugeFixingParametersQ(gfix_param)) {
        gfix_param.generic_gf.error = true;
    }

	return gfix_param;
}


/* Calculates number of sweeps to check if gauge-fixing is attained. */
static int whenNextConvergenceCheckQ(unsigned int current_sweep,
                                     const double residue, 
                                     ORGaugeFixingParameters gfix_param) {
                                 
    /*  
	 * Calls:
	 * =====
     * log10.
     * 
	 * Macros:
	 * ======
     * ESTIMATE_SWEEPS_TO_GF_PROGRESS.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * int current_sweep:   already elapsed sweeps,
	 * double residue:      index which shows how close config is to be gauge-fixed,
     * ORGaugeFixingParameters gfix_param:	overrelaxation gauge fixing parameters.
     *  
	 * Returns:
	 * =======
	 * Non-negative integer which says how long to wait to measure residue again.
     *
	 */

    if(residue < gfix_param.generic_gf.tolerance) {
        return current_sweep;
    }
    else{
        return current_sweep + 10 
                + (unsigned)(gfix_param.generic_gf.estimate_sweeps_to_gf_progress * 
                       (1.0 - log10(residue) / log10(gfix_param.generic_gf.tolerance))); 
    }

}

/* Calculates the gauge-fixing functional F, which is to be extremized when 
   performing the gauge-fixing. */
double calculateF(Mtrx3x3 * restrict U) {

    /*  
	 * Calls:
	 * =====
     * creal,
     * setNull3x3, accumulate3x3, trace3x3,
     * getLink.
     * 
     * Macros:
	 * ======
     * Nc
     * DIM, LOOP_TEMPORA, LOOP_SPATIAL, LOOP_LORENTZ_SPATIAL, T_INDX.
     * 
     * Global Variables:
     * ================
     * lattice_param
	 *
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * U:    SU(3) gluon field.
     *  
	 * Returns:
	 * =======
	 * Value of functional for the given configuration in double precision.
     *
	 */

    Mtrx3x3 U_acc;
    setNull3x3(&U_acc);

    PosIndex t;
    LOOP_TEMPORAL(t) {

        PosVec position;
        position.pos[T_INDX] = t;

        LOOP_SPATIAL(position) {
            LorentzIdx mu;
            LOOP_LORENTZ_SPATIAL(mu) {

                accumulate3x3(getLink(U, position, mu), &U_acc);

            }
        }
    }

    return creal(trace3x3(&U_acc)) / (double) (lattice_param.volume * Nc * (DIM - 1));
}

/* Calculates an alternative index to measure proximity to gauge-fixing condition.
   theta = 1/(N_c V) sum_n Re Tr[(div.A(n)).(div.A(n))^dagger]. */
double calculateTheta(Mtrx3x3 * restrict U) {

    /*  
	 * Calls:
	 * =====
     * creal,
     * hermConj3x3, prod3x3, trace3x3,
     * divergenceA.
	 *
     * Macros:
	 * ======
     * Nc,
     * LOOP_TEMPORAL, LOOP_SPATIAL, T_INDX.
     * 
     * Global Variables:
     * ================
     * lattice_param.
     * 
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * U:    SU(3) gluon field.
     *  
	 * Returns:
	 * =======
	 * Value of alternative gauge-fixing index theta in double precision.
     *
	 */

    Mtrx3x3 div_A, div_A_dagger, prod;

    double theta = 0.0;

    PosVec position;
    LOOP_TEMPORAL(position.pos[T_INDX]) {
        LOOP_SPATIAL(position) {
            divergenceA(U, position, &div_A);
            hermConj3x3(&div_A, &div_A_dagger);
            prod3x3(&div_A, &div_A_dagger, &prod);
            theta += trace3x3(&prod);
        }
    }

    return creal(theta) / (double) (Nc * lattice_param.volume);
}

/* Calculates e2, an index to measure proximity to gauge-fixing condition
   (defined in eq 6.1 of hep-lat/0301019v2). It is the normalized sum of the squares 
   of the color components of the divergence of A. */    
double calculate_e2(Mtrx3x3 * restrict U) {

    /*
	 * Calls:
	 * =====
     * decomposeAlgebraSU3,
     * divergenceA.
     * 
     * Macros:
	 * ======
     * Nc,
     * NUM_THREADS,
     * LOOP_TEMPORAL, LOOP_SPATIAL, T_INDX.
     * 
     * Global Variables:
     * ================
     * lattice_param.
	 *
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * U:    SU(3) gluon field.
     *  
	 * Returns:
	 * =======
	 * Value of gauge-fixing index e2 in double precision.
     *
	 */
    double e2 = 0.0;
    PosIndex t;
    // Parallelizing by slicing the time extent
    #pragma omp parallel for reduction (+:e2) num_threads(NUM_THREADS) schedule(dynamic) 
        LOOP_TEMPORAL(t) {
            Mtrx3x3 div_A;
            MtrxSU3Alg div_A_components;

            PosVec position;

            position.pos[T_INDX] = t;

            LOOP_SPATIAL(position) {

                divergenceA(U, position, &div_A);
                decomposeAlgebraSU3(&div_A, &div_A_components);
                for(SU3AlgIdx a = 1; a <= Nc * Nc -1; a++) {

                    /* 	Sum of the squares of the color components 
                        of the divergence of A. */
                    e2 += (div_A_components.m[a]) * (div_A_components.m[a]);

                }
            }
        }
    return e2 / (lattice_param.volume);
}

/* Accumulates the value of U_mu(n) and U_-mu(n) into w. */
static inline void accumulateFrontHearLink(Mtrx3x3 * restrict U, 
                                           const PosVec position, 
                                           LorentzIdx mu, 
                                           Mtrx3x3 * restrict w) {

    /*
	 * Calls:
	 * =====
     * conj,
     * getNeighbour,
	 * getLink.
     * 
	 * Macros:
	 * ======
     * ELEM_3X3,
     * LOOP_TEMPORAL, LOOP_SPATIAL, T_INDX.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * U: 	    SU(3) gluon field,
     * PosVec position:     specific position which is being dealt with,
     * LorentzIdx mu:       specific lorentz index which is being dealt with,
     * Mtrx3x3 * w:         matrix which holds the accumulation.
     *  
	 * Returns:
	 * =======
	 * 
	 */

    Mtrx3x3 * u_front = getLink(U,              position,             mu);
    Mtrx3x3 * u_rear  = getLink(U, getNeighbour(position, mu, REAR),  mu);

    MtrxIdx3 a, b;
    LOOP_3X3(a, b) {

        //  w += U_mu(n) + U_dagger_mu(n - mu_hat)
        
        w -> m[ELEM_3X3(a, b)] +=       u_front -> m[ELEM_3X3(a, b)]
                                 + conj(u_rear  -> m[ELEM_3X3(b, a)]);        
    
    }
}

/* Calculates w(n) = sum_mu U_mu(n) + U_dagger_mu(n-mu_hat),
   following the notation in hep-lat/9306018. */
inline static void calculate_w(Mtrx3x3 * restrict U, 
                               const PosVec position,
                               Mtrx3x3 * restrict w) {

    /*
     * Calls:
	 * =====
     * setNull3x3, accumulateFrontHearLink.
	 * 
	 * Macros:
	 * ======
     * LOOP_LORENTZ_SPATIAL.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * U: 	    SU(3) gluon field,
     * PosVec position:     specific position which is being dealt with,
     * Mtrx3x3 * w:         matrix which holds the sum of links.
     *  
	 * Returns:
	 * =======
	 * 
	 */

    setNull3x3(w);  //	Initializing w = 0

    // w calculation

    LorentzIdx mu;
    LOOP_LORENTZ_SPATIAL(mu) {

        /* w = sum_mu U_mu(n) + U_dagger_mu(n-mu_hat) */ 
        accumulateFrontHearLink(U, position, mu, w);

    }
}

/* Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
   algorithm for SU(3), with a division of the update matrix in submatrices
   following the Cabbibo-Marinari trick. Actual update is obtained after a number
   of "hits" to be performed one after another.  */
inline void updateLosAlamos(Mtrx3x3 * restrict w, 
                            const unsigned max_hits,
                            Mtrx3x3 * restrict total_update) {

    /*
	 * Calls:
	 * =====
     * setIdentity3x3, prod3x3, inverse3x3,
     * updateSubLosAlamos.
	 *
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * w:             matrix to be updated; comes from accumulating links,
     * const unsigned max_hits: iterations hits for the maximization of Tr[w(n)],
     * Submtrx sub:             index of SU(2) Cabbibo-Marinari submatrix of SU(3).
     *  
	 * Returns:
	 * =======
	 *
	 */

    Mtrx3x3 w_inv_old; 
     /* First calculates the inverse of w.
     The function will update w successively and to
     extract what was the combined update, we can 
     multiply from the right by the old inverse. */
       
    if(inverse3x3(w, &w_inv_old)) {

        /*  Local maximization is attained iteratively in SU(3),
            thus we need to make many hits ... */
        for(unsigned short hits = 1; hits <= max_hits; hits++) {

            /* ... and each hit contains the Cabbibo-Marinari subdivision. */
            for(Submtrx sub = R; sub <= T; sub++) {
                /* Submatrices are indicated by numbers from 0 to 2
                   with codenames R, S and T. */

                updateSubLosAlamos(w, sub);
                
            }
        }
    }
    else{

        /* if w has no inverse, update will be given by the identity */
        setIdentity3x3(total_update);

    }
    
    prod3x3(w, &w_inv_old, total_update);
    /*	Updates matrix with total_update. It is the
    	accumulated updates from the hits. */
}

/* Transforms gluon field U and the accumulated gauge-transformation field G by an 
   update gauge-transformation only at a given position. */
static void updateLocalUG(Mtrx3x3 * restrict U, 
                          Mtrx3x3 * restrict G, 
                          const PosVec position, 
                          const Mtrx3x3 * restrict update) {
    
    /* 
	 * Calls:
	 * =====
     * accumLeftProd3x3, accumRightProd3x3, hermConj3x3,
     * getNeighbour,
	 * getLink, getGaugetransf.
     * 
	 * Macros:
	 * ======
     * LOOP_LORENTZ.
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * U: 	    SU(3) gluon field,
     * Mtrx3x3 * G:	        SU(3) gauge-transformation field,
     * PosVec position:     specific position for the update,
     * Mtrx3x3 * update:    update gauge-transformation to be applied to the fields.
     * 
	 * Returns:
	 * =======
	 * 
	 */

    //  G'(x) = g(x) . G(x) . 1 for red-black updates

    accumLeftProd3x3(update, getGaugetransf(G, position));

    Mtrx3x3 update_dagger;
    hermConj3x3(update, &update_dagger);

    LorentzIdx mu;
    LOOP_LORENTZ(mu) {

        //	U'_mu(x) = g(x) . U_mu(x) . 1 for red-black updates

        accumLeftProd3x3(update, getLink(U, position, mu));

        //	U'_mu(x - mu) = 1 . U_mu(x - mu) . g_dagger(x) for red-black updates

        accumRightProd3x3(getLink(U, getNeighbour(position, mu, REAR), mu), 
                             &update_dagger);

    }
}

/* Performs a gauge-transformation update at the given position according to the 
   overrelaxation algorithm. Generalization of the algorithm described in 
   hep-lat/0301019v2, using the Cabbibo-Marinari submatrices trick.	*/
inline static void gaugefixOverrelaxationLocal(Mtrx3x3 * restrict U,  
                                               Mtrx3x3 * restrict G, 
                                               const PosVec position, 
                                               const ORGaugeFixingParameters params) {
    /*
     * Calls:
	 * =====
     * setIdentity3x3, power3x3Binomial,
     * calculate_w, updateLosAlamos, updateLocalUG.
	 *
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * U: 	                    SU(3) gluon field,
     * Mtrx3x3 * G:	                        SU(3) gauge-transformation field,
     * PosVec position:                     specific position for the update,
     * ORGaugeFixingParameters gfix_param:	overrelaxation gauge fixing parameters.
     * 
	 * Returns:
	 * =======
	 * 
	 */
    
    Mtrx3x3 w;
    calculate_w(U, position, &w);  

    Mtrx3x3 update_LA;
    updateLosAlamos(&w, params.hits, &update_LA);

    /*	The above function determines update_LA which would be the naïve
   	update to bring the local function to its mininum. However,	using overrelaxation,
    which means using update_LA^omega instead of  update_LA, where 1<omega<2 
    makes the converge faster. update_LA^omega is calculated using the first two terms 
    of the binomial expansion: 
    update_LA^omega=I+omega(update_LA-I)+...=I(1-omega)+omega*update_LA+... . */

    Mtrx3x3 update_OR;

    /* update_OR = update_LA^omega = Proj_SU3((I(1-omega)+omega*update_LA) */
    
    if(power3x3Binomial(&update_LA, params.omega_OR, &update_OR)){
        setIdentity3x3(&update_OR);   
        /*  If matrix could not be projected to SU3 inside
        power_3x3 binomial, then use identity as update. */
    }

    updateLocalUG(U, G, position, &update_OR);
}

/* Performs the gauge-fixing of a gluon-field following the overrelaxation 
   algorithm with the given parameters. Sets G to identity before starting the 
   procedure, overwriting anything that may have been loaded to it. */
int gaugefixOverrelaxation(Mtrx3x3 * restrict U, 
                           Mtrx3x3 * restrict G, 
                           const ORGaugeFixingParameters gfix_param) {

    /* 
	 * Calls:
	 * =====
     * fabs,
     * reunitarizeField,
     * validORGaugeFixingParametersQ, sweepsToNextMeasurement,
     * gaugefixOverrelaxationLocal, 
     * printf.
     * 
     * Macros:
	 * ======
     * LOOP_SPATIAL, 
     * BLACKRED_SWEEP, MAX_SWEEPS_TO_FIX, SWEEPS_TO_REUNITARIZATION
     * 
     * Global Variables:
     * ================
     * lattice_param
     * 
	 * Parameters:
	 * ==========
	 * Mtrx3x3 * U:                  	    SU(3) gluon field,
     * Mtrx3x3 * G:	                        SU(3) gauge-transformation field,
     * ORGaugeFixingParameters gfix_param:	overrelaxation gauge fixing parameters.
     * 
	 * Returns:
	 * =======
	 * Non-negative number of sweeps needed to gauge-fix if successful. 
     * Integer error-code otherwise. 
     *      -1: gauge-fixing algorithm didn't converge. 
     *      -2: parameters passed were invalid.
	 */

    if(!validORGaugeFixingParametersQ(gfix_param)) {
        return -2;
    }
    setFieldToIdentity(G, lattice_param.volume);

    /*	
     *  Gauge-fixing index, indicates how far we are to the Landau-gauge.
     *  It will be less than the tolerance,
     *	when the gauge-fixing is considered to be attained.
     *	Following the notation of hep-lat/0301019v2
     */
    
    //	Counter to the number of sweeps to fix config to Coulomb gauge
    int sweep = 0;
    
    /* No need to calculate residue all the time because it will typically take some 
       hundreds/thousands of sweeps to fix the gauge. */
    double new_res = gfix_param.generic_gf.gfix_proxy(U);
    printf("Sweeps: %8d.\t residue: %3.2E \n", sweep, new_res);
    
    int converge_check_due = whenNextConvergenceCheckQ(sweep, new_res, gfix_param);
    
    while(new_res > gfix_param.generic_gf.tolerance) {
        // Parallelizing by slicing the time extent
        #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic)
            for(PosIndex t = 0; t < lattice_param.n_T; t++) {
                PosVec position;
                position.pos[T_INDX] = t;
                LOOP_SPATIAL(position) {

                    BLACKRED_SWEEP(position, sweep) ?
                        /*	Implementation of the checkerboard subdivision 
                            of the lattice. */
                        
                        gaugefixOverrelaxationLocal(U, G, position, gfix_param): 0;
                        //  The actual gauge-fixing algorithm

                }
            }

        sweep++;

        if(sweep == converge_check_due) {
            double last_res = new_res;

            if(fabs((new_res = 
                    gfix_param.generic_gf.gfix_proxy(U)) - last_res)/last_res < 10E-16){
                printf("Sweeps: %8d.\t residue: %3.2E \n", sweep, new_res);
                return -1;  /* convergence too slow or not converging */
            }
            else if(new_res <= gfix_param.generic_gf.tolerance) {
                break;  /* converged */
            }
            else{
                printf("Sweeps: %8d.\t residue: %3.2E \n", sweep, new_res);
                converge_check_due = 
                                  whenNextConvergenceCheckQ(sweep, new_res, gfix_param);
            }
            
        }

        if(sweep >= gfix_param.generic_gf.max_sweeps_to_fix) {
            return -1;  /* convergence too slow or not converging */
        }

        if(!(sweep % gfix_param.generic_gf.sweeps_to_reunitarization)) {

            reunitarizeField(U, DIM * lattice_param.volume);
            reunitarizeField(G, lattice_param.volume);

        }
    }

    reunitarizeField(U, DIM * lattice_param.volume);
    reunitarizeField(G, lattice_param.volume);

    return sweep;
}