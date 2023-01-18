#ifndef GAUGEFIXING_H
#define GAUGEFIXING_H

#include <lattice.h>
#include <types.h>
#include <SU3_ops.h>

typedef struct{

    float omega_OR;
    unsigned hits;

    struct {
        double tolerance;
        double (*gfix_proxy)(Mtrx3x3 * );
        unsigned max_sweeps_to_fix;
        unsigned estimate_sweeps_to_gf_progress;
    }stop_crit;

    unsigned sweeps_to_reunitarization;

    bool error;
} ORGaugeFixingParameters;

bool validORGaugeFixingParametersQ(ORGaugeFixingParameters gfix_param);

ORGaugeFixingParameters initORGaugeFixingParameters(const char * parameter_filename);

double calculateF    (Mtrx3x3 * restrict U);
double calculateTheta(Mtrx3x3 * restrict U);
double calculate_e2  (Mtrx3x3 * restrict U);

int gaugefixOverrelaxation(Mtrx3x3 * restrict U, 
                           Mtrx3x3 * restrict G, 
                           ORGaugeFixingParameters parameters);

#endif  //GAUGEFIXING_H