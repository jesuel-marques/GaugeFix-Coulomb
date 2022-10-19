#ifndef SU3_GAUGEFIXING_PARAMETERS_H
#define SU3_GAUGEFIXING_PARAMETERS_H

//  Gauge-fixing parameters

#define TOLERANCE 1e-16              //  Tolerance for e2
#define INITIAL_SWEEPS_TO_MEASUREMENT_e2 500  //  Amount of sweeps to actually measure e2
#define SWEEPS_TO_REUNITARIZATION 250

#define MAX_HITS 2  //  Iterations hits for the maximization of Tr[w(n)]

#define OMEGA_OR 1.95  //  omega parameter of the overrelaxation

#endif