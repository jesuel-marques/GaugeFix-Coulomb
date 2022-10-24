#include <stdio.h>
#include <complex.h>
#include <omp.h>

#include  <SU3_parameters.h>

#include <lattice.h>
#include <fields.h>
#include <SU3_ops.h>

double SU3_re_tr_plaquette(mtrx_3x3 * restrict U, const pos_vec position, 
                                                  const lorentz_idx mu, 
                                                  const lorentz_idx nu){
	
    mtrx_3x3 plaquette;

    mtrx_3x3 ua, ub, uc, ud;

    pos_vec position_plus_mu = hop_pos_plus(position, mu);

    get_link_matrix(U,              position,              mu, FRONT, &ua);
    get_link_matrix(U,              position_plus_mu,      nu, FRONT, &ub);
    get_link_matrix(U, hop_pos_plus(position_plus_mu, nu), mu, REAR , &uc);
    get_link_matrix(U, hop_pos_plus(position, nu),         nu, REAR , &ud);

	prod_four_3x3(&ua, &ub, &uc, &ud, &plaquette);

    return creal(trace_3x3(&plaquette)) / Nc;

}

double spatial_plaquette_average(mtrx_3x3 * restrict U){
    //  Calculates the spatial plaquette average
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

double temporal_plaquette_average(mtrx_3x3 * restrict U){
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
                                SU3_re_tr_plaquette(U, position, T_INDX, mu);
                     
                        }    
                    }
                }
            }

            plaq_ave += plaq_ave_slice;

        }

    plaq_ave /= ((DIM - 1.0) * VOLUME);

    return plaq_ave;
}