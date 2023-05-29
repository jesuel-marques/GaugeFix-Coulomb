#include <SU3_ops.h>
#include <config_generation.h>
#include <geometry.h>

double updateLattice(Mtrx3x3* restrict U,
                     double beta,
                     double (*algorithm)(Mtrx3x3*, PosVec, LorentzIdx, double)) {
    PosIndex t;
    PosVec position;
    LorentzIdx mu;

    double deltaS = 0.0;

    LOOP_TEMPORAL(t) {
        position.pos[T_INDX] = t;
        LOOP_SPATIAL(position) {
            LOOP_LORENTZ(mu) {
                deltaS += algorithm(U, position, mu, beta);
            }
        }
    }

    return deltaS;
}