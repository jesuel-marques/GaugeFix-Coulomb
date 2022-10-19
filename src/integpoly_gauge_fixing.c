#include <tgmath.h>
#include <stdio.h>

#include  <SU3_gaugefixing_parameters.h>  //	Gauge-fixing specific parameters
#include  <SU3_parameters.h>              //	Simulation parameters

#include <SU3_ops.h>
#include <SU2_ops.h>

#include <misc.h>
#include <lattice.h>
#include <fields.h>
#include <gauge_fixing.h>
#include <matrix_power.h>
#include <integpoly_gauge_fixing.h>

#define TESTING_POWER_MATRIX
#undef TESTING_POWER_MATRIX

mtrx_3x3 average_u_temporal(mtrx_3x3 * restrict U, pos_index t){

    mtrx_3x3 u_timeslice_sum;
    mtrx_3x3 u_timeslice_ave;

    set_null_3x3(&u_timeslice_sum);


    for(pos_index k = 0; k < N_SPC; k++)
        for(pos_index j = 0; j < N_SPC; j++)
            for(pos_index i = 0; i < N_SPC; i++){
                // print_matrix_3x3(get_link(U, assign_position(i, j, k, t), T_INDX),"U",20);
                // getchar();

                accumulate_3x3(get_link(U, assign_position(i, j, k, t), T_INDX), &u_timeslice_sum);
            }
    // print_matrix_3x3(&u_timeslice_sum, "u_timeslice_sum", 16);
    scalar one_over_spatial_volume = 1.0 / pow(N_SPC, 3.0);

    // printf("%lf + I * (%lf)\n", creal(one_over_spatial_volume), cimag(one_over_spatial_volume));

    mult_by_scalar_3x3( one_over_spatial_volume, &u_timeslice_sum, &u_timeslice_ave);

    // printf("average u temporal: %d", t);
    // print_matrix_3x3(&u_timeslice_ave, "u_timeslice_ave", 16);
    return u_timeslice_ave;
}

mtrx_3x3 integ_polyakovloop(mtrx_3x3 * tempave_proj_u){
    
    mtrx_3x3 integ_polyakov_loop;
    set_identity_3x3(&integ_polyakov_loop);
    
    for(pos_index t = 0; t < N_T; t++){

        accum_left_prod_3x3(tempave_proj_u+t, &integ_polyakov_loop);

    }

    return integ_polyakov_loop;
}




int integpolyakov_gauge_fix(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G, const unsigned short config_nr) {

    printf("Integrated Polyakov Gauge-Fixing for config %d\n", config_nr);

    mtrx_3x3 u_timeslice_ave[N_T];
    mtrx_3x3 tempave_proj_u[N_T];

    mtrx_3x3 uavedag;
    mtrx_3x3 uaveproj;
    
    // OMP_PARALLEL_FOR
    for(pos_index t = 0; t < N_T; t++){


        u_timeslice_ave[t] = average_u_temporal(U, t);
        
        // printf("average u temporal: %d", t);
        // print_matrix_3x3(u_timeslice_ave+t, "", 16);
        herm_conj_3x3(u_timeslice_ave + t, &uavedag);
        SU3_CabbiboMarinari_projection(&uavedag, tempave_proj_u + t);
       
        // print_matrix_3x3(tempave_proj_u+t, "SU3 CM projected", 16);

        // getchar();
        #ifdef TESTING_POWER_MATRIX
            if(t == N_T / 2){
                print_matrix_3x3(tempave_proj_u + t, "u average middle time slice of lattice", 18);
            }
        #endif 
    }

    mtrx_3x3 P = integ_polyakovloop(tempave_proj_u);
    
    #ifdef TESTING_POWER_MATRIX
        print_matrix_3x3(&P, "P", 18);
    #endif
    mtrx_3x3 Pto1overNT, logP;
    matrix_power_3x3(&P, 1.0 / (double) N_T, &Pto1overNT);
    matrix_log_3x3(&P,  
                            &logP);
    #ifdef TESTING_POWER_MATRIX
        print_matrix_3x3(&Pto1overNT, "P to 1/Nt", 18);
        print_matrix_3x3(&logP, "log of P", 18);
    #endif

    mtrx_3x3 gt[N_T + 1 ], gdaggert[N_T + 1];
    mtrx_3x3 u_dag, aux;

    

    set_identity_3x3(gdaggert);
    set_identity_3x3(gt);
    
    set_identity_3x3(gdaggert + N_T);
    set_identity_3x3(gt + N_T);
    
    for(pos_index t = 0; t < N_T - 1 ; t++){
        
        herm_conj_3x3(tempave_proj_u + t, &u_dag); // transformar em função própria
        prod_three_3x3(&u_dag, gdaggert + t, &Pto1overNT, gdaggert + t + 1);
        herm_conj_3x3(gdaggert + t + 1, gt + t + 1);

        // #ifdef TESTING_POWER_MATRIX
        //     print_matrix_3x3(gt + t + 1, "g(t).u(t).gdag(t+1)", 16);
        // #endif
    }

    print_matrix_3x3(gt + N_T / 2, "matrix update middle", 18);

    OMP_PARALLEL_FOR
    for(pos_index t = 0; t < N_T; t++){
        mtrx_3x3 updated_u;
        pos_vec position;
        position.t = t;
        // #ifdef TESTING_POWER_MATRIX
        //     printf("t: %d %.16lf\n", t, creal(determinant_3x3(gt+t)));
        // #endif
        for(position.k = 0; position.k < N_SPC; position.k++){
            for(position.j = 0; position.j < N_SPC; position.j++){
                for(position.i = 0; position.i < N_SPC; position.i++){
                                        
                    accum_left_prod_3x3(gt + t, get_gaugetransf(G, position));
                                        
                    prod_three_3x3(gt + t, get_link(U, position, T_INDX), gdaggert + t + 1, &updated_u); //TRANSFORMAR EM FUNÇÃO
                    copy_3x3(&updated_u, get_link(U, position, T_INDX));

                    for (lorentz_idx mu = 0; mu < DIM - 1 ; mu++) {

                        prod_three_3x3(gt + t, get_link(U, position, mu), gdaggert + t, &updated_u); 
                        copy_3x3(&updated_u, get_link(U, position, mu));

                    }
                }
            }
        }
    }

}
