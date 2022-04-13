#include "../SU3_parameters.h"

#include <complex.h>
#include <omp.h>

#include "lattice.h"
#include "SU3_ops.h"


double SU3_re_tr_plaquette(double complex *U, const pos_vec position, const short mu, const short nu){
	
    double complex plaquette[3 * 3];

    double complex ua[3 * 3];
    double complex ub[3 * 3];
    double complex uc[3 * 3];
    double complex ud[3 * 3];

    pos_vec position_plus_mu = hop_position_positive(position, mu);

    get_link_matrix(U, position, mu, +1, ua);
    get_link_matrix(U, position_plus_mu, nu, +1, ub);
    get_link_matrix(U, hop_position_positive(position_plus_mu, nu), mu, -1, uc);
    get_link_matrix(U, hop_position_positive(position, nu), mu, -1, ud);

	SU3_product_four(ua, ub, uc, ud, plaquette);

    return creal(SU3_trace(plaquette));

}

double plaquette_average(double complex * U){
    double plaq_ave = 0.0;

    #pragma omp parallel for reduction (+:plaq_ave) num_threads(NUM_THREADS) schedule(dynamic) 
        // Paralelizing by slicing the time extent
        for (unsigned short t = 0; t < Nt; t++) {

            pos_vec position;

            position.t = t;
            double plaq_ave_slice = 0.0;

            for (position.i = 0; position.i < Nxyz; position.i++) {
                for (position.j = 0; position.j < Nxyz; position.j++) {
                    for (position.k = 0; position.k < Nxyz; position.k++) {
                        for (int mu = 0; mu < d; mu++){
                            for (int nu = 0; nu < d; nu++){
                                if( mu < nu){
                                
                                    plaq_ave_slice += SU3_re_tr_plaquette(U, position, mu, nu);                                    
                                
                                }
                            }
                        }
                    }
                }
            }

            plaq_ave += plaq_ave_slice;

        }

    plaq_ave /= (6.0 * Volume);

    return plaq_ave;
}