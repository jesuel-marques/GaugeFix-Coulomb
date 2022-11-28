#include <stdio.h>

#include <omp.h>

#include <fields.h>
#include <lattice.h>
#include <SU3_ops.h>
#include <SU3_parameters.h>
#include <types.h>


extern short n_SPC;
extern short n_T;

extern int volume;

double SU3_re_tr_plaquette(      Mtrx3x3 * restrict U, 
                           const PosVec position, 
                           const LorentzIdx mu, 
                           const LorentzIdx nu){
    Mtrx3x3 plaquette;

    Mtrx3x3 ua, ub, uc, ud;

    PosVec position_plus_mu = hop_pos_plus(position, mu);

    get_link_matrix(U,              position,              mu, FRONT, &ua);
    get_link_matrix(U,              position_plus_mu,      nu, FRONT, &ub);
    get_link_matrix(U, hop_pos_plus(position_plus_mu, nu), mu, REAR , &uc);
    get_link_matrix(U, hop_pos_plus(position, nu),         nu, REAR , &ud);

	prod_four_3x3(&ua, &ub, &uc, &ud, &plaquette);

    return creal(trace_3x3(&plaquette)) / Nc;

}


double spatial_plaquette_average(Mtrx3x3 * restrict U){
    //  Calculates the spatial plaquette average
    double plaq_ave = 0.0;

    #pragma omp parallel for reduction (+:plaq_ave) \
            num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (PosIndex t = 0; t < n_T; t++) {

            PosVec position;

            position.t = t;
            double plaq_ave_slice = 0.0;

            LOOP_SPATIAL(position){
                LorentzIdx mu;
                LOOP_LORENTZ_SPATIAL(mu){
                    for (LorentzIdx nu = 0; nu < DIM - 1 ; nu++){
                        if( mu < nu){

                            plaq_ave_slice += 
                                SU3_re_tr_plaquette(U, position, mu, nu);                                    
                        
                        }
                    }
                }
            }

            plaq_ave += plaq_ave_slice;

        }

    plaq_ave /= ((DIM - 1.0) * volume);

    return plaq_ave;
}


double temporal_plaquette_average(Mtrx3x3 * restrict U){
    double plaq_ave = 0.0;

    #pragma omp parallel for reduction (+:plaq_ave) \
            num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (PosIndex t = 0; t < n_T; t++) {

            PosVec position;

            position.t = t;
            double plaq_ave_slice = 0.0;

            LOOP_SPATIAL(position){
                LorentzIdx mu;
                LOOP_LORENTZ_SPATIAL(mu){
                
                    plaq_ave_slice += 
                        SU3_re_tr_plaquette(U, position, T_INDX, mu);
                
                }
            }

            plaq_ave += plaq_ave_slice;

        }

    plaq_ave /= ((DIM - 1.0) * volume);

    return plaq_ave;
}