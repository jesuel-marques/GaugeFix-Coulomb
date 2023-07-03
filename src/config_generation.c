#include <SU3_ops.h>
#include <config_generation.h>
#include <geometry.h>
#include <stdio.h>

double updateLattice(Mtrx3x3* restrict U,
                     double beta,
                     double (*algorithm)(Mtrx3x3*, PosVec, LorentzIdx, double)) {
    PosIndex t;
    PosVec position;
    LorentzIdx mu;

    double DeltaS = 0.0, deltaS = 0.0;

    int count = 0;

    LOOP_TEMPORAL(t) {
        position.pos[T_INDX] = t;
        LOOP_SPATIAL(position) {
            LOOP_LORENTZ(mu) {
                deltaS = algorithm(U, position, mu, beta);
                if (deltaS != 0.0) {
                    count++;
                }

                DeltaS += deltaS;
            }
        }
    }
    // printf("acceptance rate: %f\n", (double)count / (DIM * lattice_param.volume));
    // getchar();

    return DeltaS;
}