/*
    Definition of types used for input and output of fields.

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

#ifndef TYPES_IO_H
#define TYPES_IO_H

#include <flags.h>
#include <SU3_ops.h>
#include <tgmath.h>
#include <types.h>

/* Input scalar type */
#ifdef CONV_CFG_TO_WORKING_PRECISION
typedef complex float InScalar;
#else
typedef complex double InScalar;
#endif  //CONV_CFG_TO_WORKING_PRECISION

/* Output scalar type */
#ifdef CONV_CFG_FROM_WORKING_PRECISION
typedef complex float  OutScalar;
#else
typedef complex double OutScalar;
#endif  //CONV_CFG_FROM_WORKING_PRECISION


/* I'm aware of the rather complex structure of the input and output types below. */
/* The flags in flags.h decide which precision should be used for input and output */

#ifdef CONV_CFG_TO_WORKING_PRECISION
    /* Input configurations are in single precision */
    typedef Mtrx3x3Flt InCfgMtrx;
#else   //CONV_CFG_TO_WORKING_PRECISION
    /* Input configurations are in double precision */
    typedef Mtrx3x3Work InCfgMtrx;
#endif  //CONV_CFG_TO_WORKING_PRECISION

#ifdef CONV_GT_TO_WORKING_PRECISION
    /* Input gauge transformations are in single precision */    
    typedef Mtrx3x3Flt  InGTMtrx;
#else   //CONV_GT_TO_WORKING_PRECISION
    /* Input configurations are in double precision */
    typedef Mtrx3x3Dbl InGTMtrx;
#endif //CONV_GT_TO_WORKING_PRECISION

#ifdef CONV_CFG_FROM_WORKING_PRECISION
    /* Output configurations need to be in single precision */
    typedef Mtrx3x3Flt  OutCfgMtrx;
#else   //CONV_CFG_FROM_WORKING_PRECISION
    /* Output configurations need to be in double precision */
    typedef Mtrx3x3Dbl OutCfgMtrx;
#endif //CONV_CFG_FROM_WORKING_PRECISION

#ifdef CONV_GT_FROM_WORKING_PRECISION
    /* Output gauge-transformations need to be in single precision */
    typedef Mtrx3x3Flt  OutGTMtrx;
#else   //CONV_GT_FROM_WORKING_PRECISION
    /* Output gauge-transformations need to be in double precision */
    typedef Mtrx3x3Dbl OutGTMtrx;
#endif //CONV_GT_FROM_WORKING_PRECISION

#endif // TYPES_IO_H