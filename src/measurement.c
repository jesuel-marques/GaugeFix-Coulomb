#include <stdio.h>

#include <omp.h>

#include <fields.h>
#include <lattice.h>
#include <SU3_ops.h>

#include <types.h>

extern GeometricParameters lattice_param;

Scalar TrPlaquette(Mtrx3x3 * restrict U, 
                   const PosVec position, 
                   const LorentzIdx mu, 
                   const LorentzIdx nu) {

    /*
	 * Description:
     * ===========
	 * Calculates the trace of a specified plaquette.
     * 
	 * Calls:
	 * =====
     * creal,
     * getNeighbour, getLinkMatrix,
     * prodFour3x3, trace3x3.
	 *
	 * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * Mtrx3x3 *  U:        SU(3) gluon field,
     * PosVec position:     position to which the plaquette is attached, 
     * LorentzIdx mu:       first direction of the plaquette specification, 
     * LorentzIdx nu:       second direction of the plaquette specification.
     * 
	 * Returns:
	 * =======
     * The trace of the plaquette as a double precision complex number.
	 * 
     */

    Mtrx3x3 plaquette;

    Mtrx3x3 ua, ub, uc, ud;

    const PosVec position_plus_mu = getNeighbour(position, mu, FRONT);

    getLinkMatrix(U,              position,                     mu, FRONT, &ua);
    getLinkMatrix(U,              position_plus_mu,             nu, FRONT, &ub);
    getLinkMatrix(U, getNeighbour(position_plus_mu, nu, FRONT), mu, REAR , &uc);
    getLinkMatrix(U, getNeighbour(position,         nu, FRONT), nu, REAR , &ud);

	prodFour3x3(&ua, &ub, &uc, &ud, &plaquette);

    return trace3x3(&plaquette) / Nc;

}


Scalar averageSpatialPlaquette(Mtrx3x3 * restrict U) {

    /*
	 * Description:
     * ===========
	 * Calculates the average of the spatial plaquettes for a SU(3) field.
     * 
	 * Calls:
	 * =====
     * fprintf,
     * validGeometricParametersQ,
     * TrPlaquette.
	 *
	 * Macros:
	 * ======
     * NUM_THREADS, T_INDX, LOOP_SPATIAL, LOOP_LORENTZ_SPATIAL, DIM.
     * 
     * Global Variables:
     * ================
     * lattice_param.
     * 
	 * Parameters:
	 * ==========
     * Mtrx3x3 *  U:        SU(3) gluon field.
     * 
	 * Returns:
	 * =======
     * The average trace of the spatial plaquettes of the gluon field 
     * as a double precision complex number.
	 * 
     */

    if(!validGeometricParametersQ()) {
        fprintf(stderr, "Error in geometric parameters\n");
    }

    //  Calculates the spatial plaquette average
    Scalar plaq_ave = 0.0;
    
    // Parallelizing by slicing the time extent
    #pragma omp parallel for reduction (+:plaq_ave) \
            num_threads(NUM_THREADS) schedule(dynamic) 
        for(PosIndex t = 0; t < lattice_param.n_T; t++) {

            PosVec position;

            position.pos[T_INDX] = t;
            Scalar plaq_ave_slice = 0.0;

            LOOP_SPATIAL(position) {
                LorentzIdx mu;
                LOOP_LORENTZ_SPATIAL(mu) {
                    for(LorentzIdx nu = 0; nu < DIM - 1 ; nu++) {
                        if(mu < nu) {
                            plaq_ave_slice += 
                                TrPlaquette(U, position, mu, nu);                                    
                        }
                    }
                }
            }
            plaq_ave += plaq_ave_slice;
        }
    plaq_ave /= ((DIM - 1.0) * lattice_param.volume);

    return plaq_ave;
}


Scalar averageTemporalPlaquette(Mtrx3x3 * restrict U) {

    /*
	 * Description:
     * ===========
	 * Calculates the average of the temporal plaquettes for a SU(3) field.
     * 
	 * Calls:
	 * =====
     * fprintf,
     * validGeometricParametersQ,
     * TrPlaquette.
	 *
	 * Macros:
	 * ======
     * NUM_THREADS, LOOP_TEMPORAL, T_INDX, LOOP_SPATIAL, LOOP_LORENTZ_SPATIAL, DIM.
     * 
     * Global Variables:
     * ================
     * lattice_param.
     * 
	 * Parameters:
	 * ==========
     * Mtrx3x3 *  U:        SU(3) gluon field.
     * 
	 * Returns:
	 * =======
     * The average trace of the temporal plaquettes of the gluon field 
     * as a double precision complex number.
	 * 
     */
    
    if(!validGeometricParametersQ()) {
        fprintf(stderr, "Error in geometric parameters\n");
    }

    Scalar plaq_ave = 0.0;
    PosIndex t;
    // Parallelizing by slicing the time extent 
    #pragma omp parallel for reduction (+:plaq_ave) \
            num_threads(NUM_THREADS) schedule(dynamic) 
        LOOP_TEMPORAL(t) {

            PosVec position;

            position.pos[T_INDX] = t;
            Scalar plaq_ave_slice = 0.0;

            LOOP_SPATIAL(position) {
                LorentzIdx mu;
                LOOP_LORENTZ_SPATIAL(mu) {
                
                    plaq_ave_slice += 
                        TrPlaquette(U, position, T_INDX, mu);
                
                }
            }

            plaq_ave += plaq_ave_slice;

        }
    plaq_ave /= ((DIM - 1.0) * lattice_param.volume);

    return plaq_ave;
}