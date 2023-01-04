#ifndef GAUGEFIXING_H
#define GAUGEFIXING_H

#include <lattice.h>
#include <types.h>
#include <SU3_ops.h>

typedef struct{

    double tolerance;
    float omega_OR;
    
    bool error;
} GaugeFixingParameters;

bool validGaugeFixingParametersQ(GaugeFixingParameters gfix_param);

GaugeFixingParameters initGaugeFixingParameters(const double tolerance, 
                                                const double omega_OR);

void updateGlobalU(Mtrx3x3 * restrict U, 
                         Mtrx3x3 * restrict G);

double calculateF    (Mtrx3x3 * restrict U);
double calculateTheta(Mtrx3x3 * restrict U);
double calculate_e2  (Mtrx3x3 * restrict U);

void updateSubLosAlamos(Mtrx3x3 * restrict w, 
                              Submtrx sub); 
//GET THIS OUT OF HERE. THIS SHOULD NOT BE IN A .H FILE. ONLY HERE BECAUSE USING IN SU3_OPS

int gaugefixOverrelaxation(Mtrx3x3 * restrict U, 
                                Mtrx3x3 * restrict G, 
                                GaugeFixingParameters parameters);

#endif  //GAUGEFIXING_H