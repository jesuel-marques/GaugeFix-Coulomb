#ifndef FIELDS_H
#define FIELDS_H

#include <flags.h>
#include <lattice.h>
#include <SU3_ops.h>
#include <types.h>


//  Does the pointer arithmetic to get the correct index in the gauge transformation
#define GET_GT(G, position) G + (((position.pos[T_INDX]  * lattice_param.n_SPC \
                                 + position.pos[Z_INDX]) * lattice_param.n_SPC \
                                 + position.pos[Y_INDX]) * lattice_param.n_SPC \
                                 + position.pos[X_INDX]) 

//  Does the pointer arithmetic to get the correct index in the configuration
#define GET_LINK(U, position, mu) U + ((((position.pos[T_INDX]  * lattice_param.n_SPC \
                                        + position.pos[Z_INDX]) * lattice_param.n_SPC \
                                        + position.pos[Y_INDX]) * lattice_param.n_SPC \
                                        + position.pos[X_INDX])  * DIM                \
                                        + mu)

Mtrx3x3 *allocate3x3Field(unsigned elements);

int setFieldToIdentity(const Mtrx3x3 * restrict field,
                       unsigned elements);

int copyField(const Mtrx3x3 * restrict field, 
              unsigned elements, 
              const Mtrx3x3 * restrict field_copy);

int reunitarizeField(const Mtrx3x3 * restrict field, unsigned elements);

Scalar averageFieldDet(const Mtrx3x3 *restrict field, unsigned elements);

Mtrx3x3 * getLink(const Mtrx3x3 *U, 
                  const PosVec position,
                  const LorentzIdx mu);

void getLinkMatrix(const Mtrx3x3 * restrict U,
                   const PosVec position,
                   const LorentzIdx mu,
                   Direction dir,
                   const Mtrx3x3 * restrict u);

Mtrx3x3 *getGaugetransf(const Mtrx3x3 * restrict G,
                        const PosVec position);

void applyGaugeTransformationU(Mtrx3x3 * restrict U, 
                               Mtrx3x3 * restrict G);

#endif  //FIELDS_H