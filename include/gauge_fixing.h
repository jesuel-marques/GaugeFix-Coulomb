#ifndef GAUGEFIXING_H
#define GAUGEFIXING_H

#include <lattice.h>
#include <types.h>
#include <SU3_ops.h>

/* Generic gauge-fixing parameters */
typedef struct {

    double tolerance;
    double (*gfix_proxy)(Mtrx3x3 * );

    unsigned max_sweeps_to_fix;
    unsigned estimate_sweeps_to_gf_progress;

    unsigned sweeps_to_reunitarization;

    bool error;

} GenericGaugeFixingParameters;

/* Gauge-fixing parameters */
typedef struct{

    GenericGaugeFixingParameters generic_gf; 

    /* Overrelaxation specific gauge-fixing parameters. */
    float omega_OR;
    unsigned hits;
    
} ORGaugeFixingParameters;

bool validORGaugeFixingParametersQ(ORGaugeFixingParameters gfix_param);

ORGaugeFixingParameters initParametersORDefault();
ORGaugeFixingParameters initORGaugeFixingParameters(const char * parameter_filename);

void printORGaugeFixingParameters(ORGaugeFixingParameters gfix_param);

double calculateF    (Mtrx3x3 * restrict U);
double calculateTheta(Mtrx3x3 * restrict U);
double calculate_e2  (Mtrx3x3 * restrict U);

int gaugefixOverrelaxation(Mtrx3x3 * restrict U,
                           Mtrx3x3 * restrict G, 
                           ORGaugeFixingParameters parameters);

#endif  //GAUGEFIXING_H