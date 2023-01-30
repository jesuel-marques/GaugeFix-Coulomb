/*
    Header to lattice.c, which takes care of the geometrical aspects of lattice 
    calculations.

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

#ifndef LATTICE_H
#define LATTICE_H
#include <stdbool.h>

#include <flags.h>
#include <misc.h>
#include <types.h>

/*   Space-time lattice dimension */
#define DIM 4

/*  Associations between the numeric indices and the lorentz directions */
#define X_INDX 0
#define Y_INDX 1
#define Z_INDX 2
#define T_INDX 3


/* Test to decide if position is odd */
#define POSITION_IS_ODD(position)    (((position.pos[X_INDX]) ^ \
                                       (position.pos[Y_INDX]) ^ \
                                       (position.pos[Z_INDX]) ^ \
                                       (position.pos[T_INDX])) & 1)

/* Odd position means that the sum of the coordinates is odd and equivalente for even */

/*  if position is odd, then the XOR of the first bit of each element
    of position must be 1. Take AND with 1 select this first bit. */

/* Test to decide if position is even (not odd)*/
#define POSITION_IS_EVEN(position)   !(POSITION_IS_ODD(position))

/* Shortening for the loops in temporal coordinates */
#define LOOP_TEMPORAL(t)       for(t = 0; t < lattice_param.n_T; t++) 

/* Shortening for the loops in the spatial coordinate mu */
#define LOOP_SPATIAL_DIR(position, mu)   for(position.pos[mu] = 0; \
                                              position.pos[mu] < lattice_param.n_SPC; \
                                              position.pos[mu]++) 

/* Shortening for the loops in all spatial coordinate */                                              
#define LOOP_SPATIAL(position) LOOP_SPATIAL_DIR(position, Z_INDX)  \
                               LOOP_SPATIAL_DIR(position, Y_INDX)  \
                               LOOP_SPATIAL_DIR(position, X_INDX)

/* Shortening for the loops in temporal coordinates.  
   Parallelizing by slicing the time extent */
#define LOOP_TEMPORAL_PARALLEL(t) OMP_PARALLEL_FOR \                                
                                  LOOP_TEMPORAL(t)                           

/* Shortening for the loop in Lorentz indices */
#define LOOP_LORENTZ(mu) for(mu = 0; mu < DIM; mu++)  

/* Shortening for the loop in Lorentz indices associated to spatial directions only */
#define LOOP_LORENTZ_SPATIAL(mu) for(mu = 0; mu < DIM - 1 ; mu++)



//  Geometric types definitions

/*	Position index  */
typedef short PosIndex;

/*	Struct for position vectors */
typedef struct {
    PosIndex pos[4];
} PosVec;

/* Lorentz index */
typedef unsigned short LorentzIdx;

/*	Direction for link. A front link is the link in the positive direction from a 
    given point, whereas a rear link is a link in the negative direction  */
typedef enum {REAR, FRONT} Direction;       

/* Geometric parameters, which include the sizes of the lattice as well as the volume
   and the spatial volume for convenience. */
typedef struct {

    unsigned short n_SPC;	//   Spatial lattice size
    unsigned short n_T;		//   Temporal lattice size

    unsigned int volume;		//	Number of sites in the lattice
    unsigned int spatial_volume;	

    bool error;
} GeometricParameters;

/* lattice_param is the global lattice geometric parameters, which hold the sizes of the 
   lattice. */
extern GeometricParameters lattice_param;
/*  Although the use of a global variable here could be justified, since usually one 
    would not try to gauge-fix configurations of different sizes at the same time, it 
    is in the to-do list to actually get rid of this. */


int initGeometry(const short n_SPC, const short n_T);

bool validGeometricParametersQ(void);

bool validPositionQ(PosVec position);

bool positionmuValidQ(PosVec position, 
                       LorentzIdx mu);

PosVec assignPosition(const PosIndex x, 
                       const PosIndex y, 
                       const PosIndex z, 
                       const PosIndex t);

void printPosVec(const PosVec u);

PosVec getNeighbour(PosVec position, LorentzIdx mu, Direction dir);

#endif  //LATTICE_H