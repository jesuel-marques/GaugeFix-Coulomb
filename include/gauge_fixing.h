#ifndef GAUGEFIXING_H
#define GAUGEFIXING_H

#include <types.h>

void SU3_global_update_U(Mtrx3x3 * restrict U, 
                         Mtrx3x3 * restrict G);

double SU3_calculate_F    (Mtrx3x3 * restrict U);
double SU3_calculate_theta(Mtrx3x3 * restrict U);
double SU3_calculate_e2   (Mtrx3x3 * restrict U);

void SU3_update_sub_LosAlamos(Mtrx3x3 * restrict w, 
                              Submtrx sub); 
//GET THIS OUT OF HERE. THIS SHOULD NOT BE IN A .H FILE. ONLY HERE BECAUSE USING IN SU3_OPS

void SU3_LosAlamos(Mtrx3x3 * restrict w, 
                   Mtrx3x3 * restrict total_update);



int SU3_gaugefix_overrelaxation(Mtrx3x3 * restrict U, 
                                Mtrx3x3 * restrict G, 
                                char * config_filename, 
                                double tolerance, 
                                double omega_OR);

#endif  //GAUGEFIXING_H