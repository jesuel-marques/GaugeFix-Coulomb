#ifndef CONFIG_GENERATION_H
#define CONFIG_GENERATION_H

#include <SU3_ops.h>
#include <geometry.h>

double updateLattice(Mtrx3x3* restrict U,
                     double beta,
                     double (*algorithm)(Mtrx3x3*, PosVec, LorentzIdx, double));
#endif  // CONFIG_GENERATION_H
