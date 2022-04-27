#include "../SU3_parameters.h"

#include <stdio.h>
#include <complex.h>
#include <omp.h>

#include "lattice.h"
#include "SU3_ops.h"

double SU3_re_tr_plaquette(mtrx_3x3_double *U, const pos_vec position, 
                                               const lorentz_idx mu, 
                                               const lorentz_idx nu){
	
    mtrx_3x3_double plaquette;

    mtrx_3x3_double ua;
    mtrx_3x3_double ub;
    mtrx_3x3_double uc;
    mtrx_3x3_double ud;

    pos_vec position_plus_mu = hop_position_positive(position, mu);

    get_link_matrix(U, position,                                    mu, FRONT, &ua);
    get_link_matrix(U, position_plus_mu,                            nu, FRONT, &ub);
    get_link_matrix(U, hop_position_positive(position_plus_mu, nu), mu, REAR , &uc);
    get_link_matrix(U, hop_position_positive(position, nu),         nu, REAR , &ud);

	prod_four_3x3(&ua, &ub, &uc, &ud, &plaquette);

    return creal(trace_3x3(&plaquette))/Nc;

}

double spatial_plaquette_average(mtrx_3x3_double * U){
    double plaq_ave = 0.0;


    #pragma omp parallel for reduction (+:plaq_ave) num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {

            pos_vec position;

            position.t = t;
            double plaq_ave_slice = 0.0;

            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM - 1 ; mu++){
                            for (lorentz_idx nu = 0; nu < DIM - 1 ; nu++){
                                if( mu < nu){

                                    plaq_ave_slice += 
                                        SU3_re_tr_plaquette(U, position, mu, nu);                                    
                                
                                }
                            }
                        }
                    }
                }
            }

            plaq_ave += plaq_ave_slice;

        }

    plaq_ave /= ((DIM - 1.0) * VOLUME);

    return plaq_ave;
}

double temporal_plaquette_average(mtrx_3x3_double * U){
    double plaq_ave = 0.0;


    #pragma omp parallel for reduction (+:plaq_ave) num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (pos_index t = 0; t < N_T; t++) {

            pos_vec position;

            position.t = t;
            double plaq_ave_slice = 0.0;

            for (position.k = 0; position.k < N_SPC; position.k++) {
                for (position.j = 0; position.j < N_SPC; position.j++) {
                    for (position.i = 0; position.i < N_SPC; position.i++) {
                        for (lorentz_idx mu = 0; mu < DIM - 1 ; mu++){
                     
                                plaq_ave_slice += 
                                    SU3_re_tr_plaquette(U, position, t_index, mu);
                     
                        }    
                    }
                }
            }

            plaq_ave += plaq_ave_slice;

        }

    plaq_ave /= ((DIM - 1.0) * VOLUME);

    return plaq_ave;
}