/*
    Header to fields, which takes care of field operations on a lattice.

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

#ifndef FIELDS_H
#define FIELDS_H

#include <../settings.h>
#include <SU3_ops.h>
#include <geometry.h>
#include <types.h>

/* Does the pointer arithmetic to get the correct index for the gauge transformation */
#define GET_GT(G, position) G + (((position.pos[T_INDX] * lattice_param.n_SPC + position.pos[X_INDX]) * lattice_param.n_SPC + position.pos[Y_INDX]) * lattice_param.n_SPC + position.pos[Z_INDX])

/* Does the pointer arithmetic to get the correct index for the configuration */
#define GET_LINK(U, position, mu) U + ((((position.pos[T_INDX] * lattice_param.n_SPC + position.pos[X_INDX]) * lattice_param.n_SPC + position.pos[Y_INDX]) * lattice_param.n_SPC + position.pos[Z_INDX]) * DIM + mu)

Mtrx3x3 *allocate3x3Field(unsigned elements);

int setFieldToIdentity(const Mtrx3x3 *restrict field,
                       unsigned elements);

int copyField(const Mtrx3x3 *restrict field,
              unsigned elements,
              const Mtrx3x3 *restrict field_copy);

int reunitarizeField(const Mtrx3x3 *restrict field, unsigned elements);

Scalar averageFieldDet(const Mtrx3x3 *restrict field, unsigned elements);

Mtrx3x3 *getLink(const Mtrx3x3 *U,
                 const PosVec position,
                 const LorentzIdx mu);

void getLinkMatrix(const Mtrx3x3 *restrict U,
                   const PosVec position,
                   const LorentzIdx mu,
                   Direction dir,
                   const Mtrx3x3 *restrict u);

Mtrx3x3 *getGaugetransf(const Mtrx3x3 *restrict G,
                        const PosVec position);

void applyGaugeTransformationU(Mtrx3x3 *restrict U,
                               Mtrx3x3 *restrict G);

#endif  // FIELDS_H