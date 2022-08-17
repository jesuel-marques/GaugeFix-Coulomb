#include <tgmath.h>

#include "SU3_ops.h"
#include "SU2_ops.h"

#include "../SU3_gaugefixing_parameters.h"  //	Gauge-fixing specific parameters
#include "../SU3_parameters.h"              //	Simulation parameters

#include "lattice.h"
#include "gauge_fixing.h"
#include "integpoly_gauge_fixing.h"

mtrx_3x3 average_u_temporal(mtrx_3x3 * restrict U, pos_index t){

    mtrx_3x3 u_timeslice_sum;
    mtrx_3x3 u_timeslice_ave;

    for(pos_index k = 0; k < N_SPC; k++)
        for(pos_index j = 0; j < N_SPC; j++)
            for(pos_index i = 0; i < N_SPC; i++)
                accumulate_3x3(get_link(U, assign_position(i, j, k, t), T_INDX), &u_timeslice_sum);
    
    mult_by_scalar_3x3(1.0 / SPATIAL_VOLUME, &u_timeslice_sum, &u_timeslice_ave);

    return u_timeslice_ave;
}

work_data_type integ_polyakovloop(mtrx_3x3 * tempave_proj_u){
    
    mtrx_3x3 integ_polyakov_loop;
    set_identity_3x3(&integ_polyakov_loop);
    
    for(pos_index t = 0; t < N_T; t++){

        accum_left_prod_3x3(tempave_proj_u + t, &integ_polyakov_loop);

    }

    return trace_3x3(&integ_polyakov_loop);
}



inline static void SU3_CabbiboMarinari_projection(mtrx_3x3 * restrict u) {
    //	Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
    //	algorithm for SU(3), with a division of the update matrix in submatrices
    //	following the Cabbibo-Marinari trick. Actual update is obtained after a number
    //	of "hits" to be performed one after another.

    //  Calculates the inverse of w in the beginning.
    //  The program will update w successively and to
    //  extract what was the combined update, we can 
    //  multiply from the right by the old inverse.


    //  Local maximization is attained iteratively in SU(3),
    //  thus we need to make many hits ...
    for (unsigned short hits = 1; hits <= 8; hits++) {

        //	... and each hit contains the Cabbibo-Marinari subdivision
        for (submatrix sub = R; sub <= T; sub++) {
            //	Submatrices are indicated by numbers from 0 to 2
            //  with codenames R, S and T

            SU3_update_sub_LosAlamos(u, sub);
            
        }
    }
}

int integpolyakov_gauge_fix(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G, const unsigned short config_nr) {

    printf("Integrated Polyakov Gauge-Fixing for config %d\n", config_nr);

    mtrx_3x3 u_timeslice_ave[N_T];
    mtrx_3x3 tempave_proj_u[N_T];
    
    omp_parallel_for
    for(pos_index t = 0; t < N_T; t++){


        u_timeslice_ave[t] = average_u_temporal(U, t);
        
        printf("average u temporal: %d", t);
        print_matrix_3x3(u_timeslice_ave+t, "", 16);

        SU3_CabbiboMarinari_projection(u_timeslice_ave + t);
        print_matrix_3x3(u_timeslice_ave+t, "SU3 CM projected", 16);

        getchar();
    }



    double Pto1overNT = pow(integ_polyakovloop(u_timeslice_ave), 1.0 / N_T);

    mtrx_3x3 gt[N_T], gdaggert[N_T];
    mtrx_3x3 u_dag, aux;

    set_identity_3x3(gdaggert);
    set_identity_3x3(gt);

    for(pos_index t = 0; t < N_T - 1 ; t++){
        
        herm_conj_3x3(u_timeslice_ave + t, &u_dag); // transformar em função própria
        prod_3x3(&u_dag, gdaggert + t, gdaggert + t + 1);
        mult_by_scalar_3x3(Pto1overNT, gdaggert + t + 1, &aux);  // fazer multiplicação com acumulação
        copy_3x3(&aux, gdaggert + t + 1);
        herm_conj_3x3(gdaggert + t + 1, gt + t + 1);
    }
    omp_parallel_for
    for(pos_index t = 0; t < N_T; t++){
        pos_vec position;
        position.t = t;
        for(position.k = 0; position.k < N_SPC; position.k++){
            for(position.j = 0; position.j < N_SPC; position.j++){
                for(position.i = 0; position.i < N_SPC; position.i++){

                    accum_left_prod_3x3(gt + t, get_gaugetransf(G,position));

                }
            }
        }
    }

}
