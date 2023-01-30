/*
    Header to gauge_fixing.c, which contain routines to gauge fix a configuration.

    Copyright (C) 2023  Jesuel Marques

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Contact: jesuel.leal@usp.br
    
 */

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