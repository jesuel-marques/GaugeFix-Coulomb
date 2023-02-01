/*
    Sets flags related to format conversions of input/output fields.

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

#ifndef FLAGS_H
#define FLAGS_H

/* Uncomment if conversion from big to little endian is needed in config being read */
#define     NEED_BYTE_SWAP_IN                

/* Uncomment if conversion from big to little endian is needed in config being written*/
// #define  NEED_BYTE_SWAP_OUT                  


/* Uncomment if reading config in float and need to convert to double precision */
// #define     CONV_CFG_TO_WORKING_PRECISION

/*  Uncomment if writing config in float and need to convert from double precision */
// #define     CONV_CFG_FROM_WORKING_PRECISION


/*  Uncomment if reading gauge-transf in float and 
    need to convert to double precision */
// #define  CONV_GT_TO_WORKING_PRECISION

/*  Uncomment if writing gauge-transf in float and 
    need to convert from double precision */
// #define  CONV_GT_FROM_WORKING_PRECISION


// #define CHECK_POSITION_BOUNDS
//  Uncomment to check that only positions within bounds are being accessed

#endif  //FLAGS_H