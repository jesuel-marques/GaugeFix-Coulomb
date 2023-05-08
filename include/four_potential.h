/*
    header to four_potential.c, which contains routines to calculate the gauge
    four-potential field and related quantities.

    Copyright (C) 2023  Jesuel Marques

    This program is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this
    program.  If not, see <https://www.gnu.org/licenses/>.

    Contact: jesuel.leal@usp.br

 */

#ifndef FOURVECTORFIELD_H
#define FOURVECTORFIELD_H

#include <geometry.h>
#include <SU3_ops.h>

typedef enum {SPATIAL, QUADRI} DivergenceType;

void calculateA (Mtrx3x3 * restrict U, 
                 const PosVec position,
                 const LorentzIdx mu, 
                 Mtrx3x3 * restrict A);

void DivergenceA(Mtrx3x3 * restrict U,
                 const PosVec position, 
                 Mtrx3x3 * restrict div_A,
                 DivergenceType divergence_type);

#endif  //FOURVECTORFIELD_H