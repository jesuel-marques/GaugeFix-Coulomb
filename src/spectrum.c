#include <SU3_ops.h>
#include <dirac.h>
#include <geometry.h>
#include <string.h>

void calculateCorrelator(Scalar* inverse[4][Nc], double* correlator, const char* type) {
    PosVec position;

    DiracIdx alpha, beta;
    MtrxIdx3 a, b;

    if (!strcmp(type, "pion")) {
        LOOP_TEMPORAL(position.pos[T_INDX]) {
            *correlator = 0.0;
            LOOP_SPATIAL(position) {
                LOOP_DIRAC(alpha) {
                    LOOP_DIRAC(beta) {
                        LOOP_3(a) {
                            LOOP_3(b) {
                                *correlator += pow(cabs(*(ELEM_VEC_POSDC(inverse[beta][b], position, alpha, a))), 2.0);
                            }
                        }
                    }
                }
            }
            correlator++;
        }
    }
}