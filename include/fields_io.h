/*
    Header to fields_io.c, which takes care of input and output of lattice fields.

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

#ifndef FIELDS_IO_H
#define FIELDS_IO_H

#include <../settings.h>
#include <types_io.h>
#include <geometry.h>

InGTMtrx  *getGaugetransfIn (InGTMtrx  * restrict G_in,  
                               const PosVec position);
OutGTMtrx *getGaugetransfOut(OutGTMtrx * restrict G_out, 
                               const PosVec position);

InCfgMtrx  *getLinkIn (InCfgMtrx  * restrict U_in,
                         const PosVec position, 
                         const LorentzIdx mu);
OutCfgMtrx *getLinkOut(OutCfgMtrx * restrict U_out,
                         const PosVec position,
                         const LorentzIdx mu);

int loadConfig (const Mtrx3x3 * restrict U,
                     char * config_filename),
    writeConfig(Mtrx3x3 * restrict U,
                     char * config_filename);

int loadGaugeTransf (Mtrx3x3 * restrict G,
                           char * gauge_transf_filename),
    writeGaugeTransf(Mtrx3x3 * restrict G,
                           char * gauge_transf_filename);

#endif  //FIELDS_IO_H