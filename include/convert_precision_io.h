/*
    Header to convert_precision_io.c, which takes care of converting the input or output
    fields to the format other programs in this suite use internally.

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

#ifndef CONVERT_PRECISION_IO_H
#define  CONVERT_PRECISION_IO_H


#include <../settings.h>
#include <types_io.h>
#include <geometry.h>
#include <types.h>
#include <SU3_ops.h>


#ifdef CONV_CFG_TO_WORKING_PRECISION

void convertCfg_in_work(InCfgMtrx * restrict U_in, Mtrx3x3Work  * restrict U_work);

#endif  //CONV_CFG_TO_WORKING_PRECISION

#ifdef CONV_GT_TO_WORKING_PRECISION

void convertGT_in_work(InGTMtrx * restrict G_in, Mtrx3x3Work * restrict G_work);

#endif  //CONV_GT_TO_WORKING_PRECISION


#ifdef CONV_CFG_FROM_WORKING_PRECISION

void convertCfg_work_out(Mtrx3x3Work * restrict U_work, OutCfgMtrx * restrict U_out);

#endif  //CONV_CFG_FROM_WORKING_PRECISION

#ifdef CONV_GT_FROM_WORKING_PRECISION

void convertGT_work_out(Mtrx3x3Work *G_work, OutGTMtrx *G_out);

#endif  //CONV_GT_FROM_WORKING_PRECISION

#endif  //CONVERT_PRECISION_IO_H