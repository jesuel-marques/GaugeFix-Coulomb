/*
    Definitions of types used internally in the code.

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

#ifndef TYPES_H
#define TYPES_H

#include <tgmath.h>

//  Data types definitions

typedef complex double WorkScalarType;      //  Calculations internally done in double
typedef WorkScalarType Scalar;              //  Default scalar is a scalar of type Work

#endif  //TYPES_H