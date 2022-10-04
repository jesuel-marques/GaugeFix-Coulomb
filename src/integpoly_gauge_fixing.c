#include <tgmath.h>
#include <stdio.h>

#include  <SU3_gaugefixing_parameters.h>  //	Gauge-fixing specific parameters
#include  <SU3_parameters.h>              //	Simulation parameters

#include <SU3_ops.h>
#include <SU2_ops.h>

#include <lattice.h>
#include <gauge_fixing.h>
#include <matrix_power.h>
#include <integpoly_gauge_fixing.h>


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
    work_data_type one_over_spatial_volume = 1.0 / pow(N_SPC, 3.0);

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



static void SU3_CabbiboMarinari_projection(mtrx_3x3 * restrict w, 
                                              mtrx_3x3 * restrict total_update) {
    //	Calculates the update matrix A from w(n)=g(n).h(n) as in the Los Alamos
    //	algorithm for SU(3), with a division of the update matrix in submatrices
    //	following the Cabbibo-Marinari trick. Actual update is obtained after a number
    //	of "hits" to be performed one after another.

    mtrx_3x3 w_inv_old; 
    //  Calculates the inverse of w in the beginning.
    //  The program will update w successively and to
    //  extract what was the combined update, we can 
    //  multiply from the right by the old inverse.

       
    if(inverse_3x3(w, &w_inv_old)){

        //  Local maximization is attained iteratively in SU(3),
        //  thus we need to make many hits ...
        for (unsigned short hits = 1; hits <= 300; hits++) {

            //	... and each hit contains the Cabbibo-Marinari subdivision
            for (submatrix sub = R; sub <= T; sub++) {
                //	Submatrices are indicated by numbers from 0 to 2
                //  with codenames R, S and T

                SU3_update_sub_LosAlamos(w, sub);
                // projection_SU3(w);
                // printf("re Tr w: %.20lf\n", creal(trace_3x3(w)));
                
            }
        }
    }
    else{
        //  if w has no inverse, update will be given by the identity
        set_identity_3x3(total_update);
    }
    
    prod_3x3(w, &w_inv_old, total_update);
    //	Updates matrix to total_update. It is the
    //	accumulated updates from the hits.
}

int integpolyakov_gauge_fix(mtrx_3x3 * restrict U, mtrx_3x3 * restrict G, const unsigned short config_nr) {

    printf("Integrated Polyakov Gauge-Fixing for config %d\n", config_nr);

    mtrx_3x3 u_timeslice_ave[N_T];
    mtrx_3x3 tempave_proj_u[N_T];

    mtrx_3x3 uavedag;
    mtrx_3x3 uaveproj;
    
    // omp_parallel_for
    for(pos_index t = 0; t < N_T; t++){


        u_timeslice_ave[t] = average_u_temporal(U, t);
        
        // printf("average u temporal: %d", t);
        // print_matrix_3x3(u_timeslice_ave+t, "", 16);
        herm_conj_3x3(u_timeslice_ave + t, &uavedag);
        SU3_CabbiboMarinari_projection(&uavedag, tempave_proj_u + t);
        // print_matrix_3x3(tempave_proj_u+t, "SU3 CM projected", 16);

        // getchar();
    }

    mtrx_3x3 P = integ_polyakovloop(tempave_proj_u);
    print_matrix_3x3(&P, "P", 16);
    mtrx_3x3 Pto1overNT;
    matrix_power_3x3(&P, 1.0 / (double) N_T, &Pto1overNT);
    print_matrix_3x3(&P, "P to 1/Nt", 16);

    mtrx_3x3 gt[N_T], gdaggert[N_T];
    mtrx_3x3 u_dag, aux;

    set_identity_3x3(gdaggert);
    set_identity_3x3(gt);
    mtrx_3x3 result;
    for(pos_index t = 0; t < N_T - 1 ; t++){
        
        herm_conj_3x3(tempave_proj_u + t, &u_dag); // transformar em função própria
        prod_three_3x3(&u_dag, gdaggert + t, &Pto1overNT, gdaggert + t + 1);
        herm_conj_3x3(gdaggert + t + 1, gt + t + 1);

        prod_three_3x3(gt, tempave_proj_u, gdaggert + 1, &result);
        print_matrix_3x3(&result, "g(t).u(t).gdag(t+1)", 16);
        getchar();
    }

    omp_parallel_for
    for(pos_index t = 0; t < N_T; t++){
        pos_vec position;
        position.t = t;
        printf("t: %d %.16lf\n", t, creal(determinant_3x3(gt+t)));
        for(position.k = 0; position.k < N_SPC; position.k++){
            for(position.j = 0; position.j < N_SPC; position.j++){
                for(position.i = 0; position.i < N_SPC; position.i++){
                    
                    accum_left_prod_3x3(gt + t, get_gaugetransf(G, position));

                }
            }
        }
    }

    SU3_global_update_U(U, G);

}