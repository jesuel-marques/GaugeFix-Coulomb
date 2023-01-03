#ifndef GAUGEFIXING_H
#define GAUGEFIXING_H

#include <lattice.h>
#include <types.h>
#include <SU3_ops.h>

typedef struct{

    double tolerance;
    float omega_OR;
    
    bool error;
} gauge_fixing_parameters;


gauge_fixing_parameters init_gaugefixing_parameters(const double tolerance, 
                                                    const double omega_OR);

void SU3_global_update_U(Mtrx3x3 * restrict U, 
                         Mtrx3x3 * restrict G);

double SU3_calculate_F    (Mtrx3x3 * restrict U);
double SU3_calculate_theta(Mtrx3x3 * restrict U);
double SU3_calculate_e2   (Mtrx3x3 * restrict U);

void SU3_update_sub_LosAlamos(Mtrx3x3 * restrict w, 
                              Submtrx sub); 
//GET THIS OUT OF HERE. THIS SHOULD NOT BE IN A .H FILE. ONLY HERE BECAUSE USING IN SU3_OPS

int SU3_gaugefix_overrelaxation(Mtrx3x3 * restrict U, 
                                Mtrx3x3 * restrict G, 
                                gauge_fixing_parameters parameters);

#endif  //GAUGEFIXING_H