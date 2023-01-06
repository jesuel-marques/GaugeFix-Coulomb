#include <stdio.h>

#include <omp.h>

#include <fields.h>
#include <lattice.h>
#include <SU3_ops.h>

#include <types.h>

extern GeometricParameters lattice_param;

double ReTrPlaquette(Mtrx3x3 * restrict U, 
                     const PosVec position, 
                     const LorentzIdx mu, 
                     const LorentzIdx nu) {
    Mtrx3x3 plaquette;

    Mtrx3x3 ua, ub, uc, ud;

    const PosVec position_plus_mu = getNeighbour(position, mu, FRONT);

    getLinkMatrix(U,              position,                     mu, FRONT, &ua);
    getLinkMatrix(U,              position_plus_mu,             nu, FRONT, &ub);
    getLinkMatrix(U, getNeighbour(position_plus_mu, nu, FRONT), mu, REAR , &uc);
    getLinkMatrix(U, getNeighbour(position,         nu, FRONT), nu, REAR , &ud);

	prodFour3x3(&ua, &ub, &uc, &ud, &plaquette);

    return creal(trace3x3(&plaquette)) / Nc;

}


double averageSpatialPlaquette(Mtrx3x3 * restrict U) {

    if(!validGeometricParametersQ()){
        fprintf(stderr, "Error in geometric parameters\n");
    }

    //  Calculates the spatial plaquette average
    double plaq_ave = 0.0;
    
    // Parallelizing by slicing the time extent
    #pragma omp parallel for reduction (+:plaq_ave) \
            num_threads(NUM_THREADS) schedule(dynamic) 
        for(PosIndex t = 0; t < lattice_param.n_T; t++) {

            PosVec position;

            position.pos[T_INDX] = t;
            double plaq_ave_slice = 0.0;

            LOOP_SPATIAL(position) {
                LorentzIdx mu;
                LOOP_LORENTZ_SPATIAL(mu) {
                    for(LorentzIdx nu = 0; nu < DIM - 1 ; nu++) {
                        if( mu < nu) {

                            plaq_ave_slice += 
                                ReTrPlaquette(U, position, mu, nu);                                    
                        
                        }
                    }
                }
            }

            plaq_ave += plaq_ave_slice;

        }

    plaq_ave /= ((DIM - 1.0) * lattice_param.volume);

    return plaq_ave;
}


double averageTemporalPlaquette(Mtrx3x3 * restrict U) {
    
    if(!validGeometricParametersQ()){
        fprintf(stderr, "Error in geometric parameters\n");
    }

    double plaq_ave = 0.0;
    PosIndex t;
    // Parallelizing by slicing the time extent 
    #pragma omp parallel for reduction (+:plaq_ave) \
            num_threads(NUM_THREADS) schedule(dynamic) 
        LOOP_TEMPORAL(t) {

            PosVec position;

            position.pos[T_INDX] = t;
            double plaq_ave_slice = 0.0;

            LOOP_SPATIAL(position) {
                LorentzIdx mu;
                LOOP_LORENTZ_SPATIAL(mu) {
                
                    plaq_ave_slice += 
                        ReTrPlaquette(U, position, T_INDX, mu);
                
                }
            }

            plaq_ave += plaq_ave_slice;

        }

    plaq_ave /= ((DIM - 1.0) * lattice_param.volume);

    return plaq_ave;
}