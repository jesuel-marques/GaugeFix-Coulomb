#include <stdio.h>
#include <tgmath.h>

#include <integpoly_gauge_fixing.h>
#include <lattice.h>
#include <matrix_power.h>
#include <misc.h>
#include <SU3_ops.h>
#include <types.h>


#define TESTING_POWER_MATRIX
#undef TESTING_POWER_MATRIX

extern short n_SPC;
extern short n_T;

Mtrx3x3 average_u_temporal(Mtrx3x3 * restrict U, 
                           PosIndex t){
    Mtrx3x3 u_timeslice_sum;
    Mtrx3x3 u_timeslice_ave;

    set_null_3x3(&u_timeslice_sum);

    PosVec position = {.t = t};
    LOOP_SPATIAL(position){

        accumulate_3x3(get_link(U, position, T_INDX), &u_timeslice_sum);
    }
    Scalar one_over_spatial_volume = 1.0 / pow(n_SPC, 3.0);


    mult_by_scalar_3x3( one_over_spatial_volume, &u_timeslice_sum, &u_timeslice_ave);

    return u_timeslice_ave;
}


Mtrx3x3 integ_polyakovloop(Mtrx3x3 * tempave_proj_u){
    Mtrx3x3 integ_polyakov_loop;
    set_identity_3x3(&integ_polyakov_loop);
    
    for(PosIndex t = 0; t < n_T; t++){

        accum_left_prod_3x3(tempave_proj_u+t, &integ_polyakov_loop);

    }

    return integ_polyakov_loop;
}


int integ_polyakov_gaugefix(      Mtrx3x3 * restrict U, 
                                  Mtrx3x3 * restrict G) {

    Mtrx3x3 u_timeslice_ave[n_T];
    Mtrx3x3 tempave_proj_u[n_T];

    Mtrx3x3 uavedag;
    Mtrx3x3 uaveproj;
    
    // OMP_PARALLEL_FOR
    for(PosIndex t = 0; t < n_T; t++){


        u_timeslice_ave[t] = average_u_temporal(U, t);
        
        herm_conj_3x3(u_timeslice_ave + t, &uavedag);
        SU3_CabbiboMarinari_projection(&uavedag, tempave_proj_u + t);
     
        return 0;
    }

    Mtrx3x3 P = integ_polyakovloop(tempave_proj_u);
    
    Mtrx3x3 Pto1overNT, logP;
    matrix_power_3x3(&P, 1.0 / (double) n_T, &Pto1overNT);
    matrix_log_3x3(&P,  
                            &logP);
    
    Mtrx3x3 gt[n_T + 1], gdaggert[n_T + 1];
    Mtrx3x3 u_dag, aux;

    set_identity_3x3(gdaggert);
    set_identity_3x3(gt);
    
    set_identity_3x3(gdaggert + n_T);
    set_identity_3x3(gt + n_T);
    
    for(PosIndex t = 0; t < n_T - 1 ; t++){
        
        herm_conj_3x3(tempave_proj_u + t, &u_dag); // transformar em função própria
        prod_three_3x3(&u_dag, gdaggert + t, &Pto1overNT, gdaggert + t + 1);
        herm_conj_3x3(gdaggert + t + 1, gt + t + 1);

        
    }


    OMP_PARALLEL_FOR
    for(PosIndex t = 0; t < n_T; t++){
        Mtrx3x3 updated_u;
        PosVec position;
        position.t = t;
        
        LOOP_SPATIAL(position){
                                        
            accum_left_prod_3x3(gt + t, get_gaugetransf(G, position));
                                
            prod_three_3x3(gt + t, get_link(U, position, T_INDX), gdaggert + t + 1, 
                           &updated_u); //TRANSFORMAR EM FUNÇÃO
            copy_3x3(&updated_u, get_link(U, position, T_INDX));
            LorentzIdx mu;

            LOOP_LORENTZ_SPATIAL(mu){

                prod_three_3x3(gt + t, get_link(U, position, mu), gdaggert + t, 
                               &updated_u); 
                copy_3x3(&updated_u, get_link(U, position, mu));

            }
        }
    }

    return 0;
}