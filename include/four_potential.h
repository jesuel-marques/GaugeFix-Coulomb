#ifndef FOURVECTORFIELD_H
#define FOURVECTORFIELD_H

#include <lattice.h>
#include <SU3_ops.h>

void calculateA (Mtrx3x3 * restrict U, 
                 const PosVec position,
                 const LorentzIdx mu, 
                 Mtrx3x3 * restrict A);

void divergenceA(Mtrx3x3 * restrict U,
                 const PosVec position, 
                 Mtrx3x3 * restrict div_A);

#endif  //FOURVECTORFIELD_H