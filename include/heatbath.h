/*
    Header to heatbath.c, which contains functions used in the heatbath algorithm for
    generating gauge configurations.

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

#ifndef HEATBATH_H
#define HEATBATH_H

#include <SU3_ops.h>
#include <geometry.h>

double MetropolisSU3(Mtrx3x3* restrict U, PosVec position, LorentzIdx mu, double beta);
double HeatBathSU3(Mtrx3x3* U, PosVec position, LorentzIdx mu, double beta);

#endif  // HEATBATH_H