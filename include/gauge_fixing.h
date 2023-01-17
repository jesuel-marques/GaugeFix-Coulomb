#ifndef GAUGEFIXING_H
#define GAUGEFIXING_H

#include <lattice.h>
#include <types.h>
#include <SU3_ops.h>

typedef struct{

    double tolerance;
    float omega_OR;
    
    bool error;
} ORGaugeFixingParameters;

bool validORGaugeFixingParametersQ(ORGaugeFixingParameters gfix_param);

ORGaugeFixingParameters initORGaugeFixingParameters(const double tolerance, 
                                                    const double omega_OR);


double calculateF    (Mtrx3x3 * restrict U);
double calculateTheta(Mtrx3x3 * restrict U);
double calculate_e2  (Mtrx3x3 * restrict U);

int gaugefixOverrelaxation(Mtrx3x3 * restrict U, 
                           Mtrx3x3 * restrict G, 
                           ORGaugeFixingParameters parameters);

#endif  //GAUGEFIXING_H