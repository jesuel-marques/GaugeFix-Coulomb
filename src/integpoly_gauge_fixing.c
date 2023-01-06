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

extern GeometricParameters lattice_param;

Mtrx3x3 average_u_Temporal(Mtrx3x3 * restrict U, 
                           PosIndex t) {
    Mtrx3x3 u_timeslice_sum;
    Mtrx3x3 u_timeslice_ave;

    setNull3x3(&u_timeslice_sum);

    PosVec position = {.pos={t,0,0,0}};
    LOOP_SPATIAL(position) {

        accumulate3x3(getLink(U, position, T_INDX), &u_timeslice_sum);
    }
    Scalar one_over_spatial_volume = 1.0 / pow(lattice_param.n_SPC, 3.0);


    multByScalar3x3( one_over_spatial_volume, &u_timeslice_sum, &u_timeslice_ave);

    return u_timeslice_ave;
}


Mtrx3x3 integPolyakovLoop(Mtrx3x3 * tempave_proj_u) {
    Mtrx3x3 integ_polyakov_loop;
    setIdentity3x3(&integ_polyakov_loop);
    
    for(PosIndex t = 0; t < lattice_param.n_T; t++) {

        accumLeftProd3x3(tempave_proj_u+t, &integ_polyakov_loop);

    }

    return integ_polyakov_loop;
}


int integPolyakovGaugefix(Mtrx3x3 * restrict U, 
                          Mtrx3x3 * restrict G) {

    Mtrx3x3 u_timeslice_ave[lattice_param.n_T];
    Mtrx3x3 tempave_proj_u[lattice_param.n_T];

    Mtrx3x3 uavedag;
    Mtrx3x3 uaveproj;
    
    // OMP_PARALLEL_FOR
    for(PosIndex t = 0; t < lattice_param.n_T; t++) {


        u_timeslice_ave[t] = average_u_Temporal(U, t);
        
        hermConj3x3(u_timeslice_ave + t, &uavedag);
        projectSU3CabbiboMarinari(&uavedag, tempave_proj_u + t);
     
        return 0;
    }

    Mtrx3x3 P = integPolyakovLoop(tempave_proj_u);
    
    Mtrx3x3 Pto1overNT, logP;
    powerMtrx3x3(&P, 1.0 / (double) lattice_param.n_T, &Pto1overNT);
    logMtrx3x3(&P,  
                            &logP);
    
    Mtrx3x3 gt[lattice_param.n_T + 1], gdaggert[lattice_param.n_T + 1];
    Mtrx3x3 u_dag, aux;

    setIdentity3x3(gdaggert);
    setIdentity3x3(gt);
    
    setIdentity3x3(gdaggert + lattice_param.n_T);
    setIdentity3x3(gt + lattice_param.n_T);
    
    for(PosIndex t = 0; t < lattice_param.n_T - 1 ; t++) {
        
        hermConj3x3(tempave_proj_u + t, &u_dag); // transformar em função própria
        prodThree3x3(&u_dag, gdaggert + t, &Pto1overNT, gdaggert + t + 1);
        hermConj3x3(gdaggert + t + 1, gt + t + 1);

        
    }


    OMP_PARALLEL_FOR
    for(PosIndex t = 0; t < lattice_param.n_T; t++) {
        Mtrx3x3 updated_u;
        PosVec position;
        position.pos[T_INDX] = t;
        
        LOOP_SPATIAL(position) {
                                        
            accumLeftProd3x3(gt + t, getGaugetransf(G, position));
                                
            prodThree3x3(gt + t, getLink(U, position, T_INDX), gdaggert + t + 1, 
                           &updated_u); //TRANSFORMAR EM FUNÇÃO
            copy3x3(&updated_u, getLink(U, position, T_INDX));
            LorentzIdx mu;

            LOOP_LORENTZ_SPATIAL(mu) {

                prodThree3x3(gt + t, getLink(U, position, mu), gdaggert + t, 
                               &updated_u); 
                copy3x3(&updated_u, getLink(U, position, mu));

            }
        }
    }

    return 0;
}