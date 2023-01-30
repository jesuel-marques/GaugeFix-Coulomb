/*
    lattice takes care of the geometrical aspects of lattice calculations.

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

#include <stdio.h>
#include <stdbool.h>

#include <flags.h>
#include <lattice.h>
#include <SU3_ops.h>
#include <types.h>

GeometricParameters lattice_param;

/* Modulo function which deals with negative numbers. */
static inline int mod(short a, unsigned short b) {
    
    /*  
	 * Calls:
	 * =====
     *
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
     * short a:     number whose mod is being taken,
     * short b:     the mod base.
     * 
	 * Returns:
	 * =======
	 * a mod b. 
     * 
     */

    short r = a % b;
    return r >= 0 ? r : r + b;
}

/* Implements the periodic boundary conditions at the geometric level. */
PosVec makePeriodicBound(PosVec position) {

    /*
	 * Calls:
	 * =====
     * mod.
     * 
     * Macros:
	 * ======
     * LOOP_LORENTZ, T_INDX.
     * 
     * Global Variables:
     * ================
     * lattice_param
     *  
	 * Parameters:
	 * ==========
     * PosVec position:    position to be transformed into a valid lattice size in 
     *                     the correct range and taking into account the periodic 
     *                     boundary conditions.
     * 
	 * Returns:
	 * =======
	 * position with fields in the range between 0 and N_\mu-1 for all directions.
     * 
     */

    PosVec valid_position;
    
    // if(!validGeometricParametersQ()) {
    //     fprintf(stderr, "Error in geometric parameters\n";
    //     return position;
    // }
    LorentzIdx mu;
    LOOP_LORENTZ(mu) {
        valid_position.pos[mu] = mod(position.pos[mu], mu != T_INDX ? 
                                                       lattice_param.n_SPC :
                                                       lattice_param.n_T);
    }
    
    return valid_position;   
}


/* Returns the front or rear (depending on dir) neighbour of position in Lorentz 
   direction mu.  */
inline PosVec getNeighbour(PosVec position, LorentzIdx mu, Direction dir) {

    /* 
	 * Calls:
	 * =====
     * printf,
     * makePeriodicBound.
     * 
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     * 
	 * Parameters:
	 * ==========
     * PosVec position:     position whose neighbour is desired,
     * LorentzIdx mu:       Lorentz direction of the neighbour, 
     * Direction dir:       variable to decide if one wants the front or rear neighbour.
     * 
	 * Returns:
	 * =======
     * A struct PosVec with the neighbour position.
     */

    PosVec neighbour = position; 

    if(dir == REAR) {
        neighbour.pos[mu]--;
    }
    else{
        neighbour.pos[mu]++;
    } 

    return makePeriodicBound(neighbour); 
    
}


/* Check if current geometric parameters are valid. */
bool validGeometricParametersQ(void) {

    /*
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * DIM
     * 
     * Global Variables:
     * ================
     * lattice_param
     *  
	 * Parameters:
	 * ==========
     * 
	 * Returns:
	 * =======
	 * 
     */

    short ns = lattice_param.n_SPC;
    short nt = lattice_param.n_T;

    return (lattice_param.n_SPC          >  0                && 
            lattice_param.n_T            >  0                &&
            lattice_param.error          != 1                &&
            lattice_param.spatial_volume == ns * ns * ns     &&
            lattice_param.volume         == ns * ns * ns * nt  );
}


/* Initializes the global variable lattice_param which contains the lattice 
   geometric parameters, performing check to see if the parameters set are
   reasonable. After that it generates a table of neighbours for neighbour lookups. */
int initGeometry(const short n_s, const short n_t) {

    /*  
	 * Calls:
	 * =====
     * validGeometricParametersQ.
     * 
     * Macros:
	 * ======
     * DIM
     * 
     * Global Variables:
     * ================
     * lattice_param
     *  
	 * Parameters:
	 * ==========
     * short n_s:   lattice spatial extent
     * short n_t:   lattice temporal extent
     * 
	 * Returns:
	 * =======
	 * 0 on success, i. e. if the geometric parameters passed were valid, 1 on failure.
     * 
     */

    lattice_param.n_SPC 		= n_s;
	lattice_param.n_T   		= n_t;

	if(lattice_param.n_SPC <= 0 || lattice_param.n_T <= 0) {
		lattice_param.error = 1;
	}

	lattice_param.spatial_volume = lattice_param.n_SPC * 
                                   lattice_param.n_SPC * 
                                   lattice_param.n_SPC;	
	lattice_param.volume 	     = lattice_param.spatial_volume * 
                                   lattice_param.n_T;

    if(!validGeometricParametersQ()) {
        return 1;
    }

    return 0;
}


/* Checks if all position components are in the allowed range between 0 and N_mu. */
inline bool validPositionQ(PosVec position) {

    /*  
	 * Calls:
	 * =====
     * 
     * Macros:
	 * ======
     * LOOP_LORENTZ, T_INDX.
     * 
     * Global Variables:
     * ================
     * lattice_param
     *  
	 * Parameters:
	 * ==========
     * PosVec position:     position to be checked.
     * 
	 * Returns:
	 * =======
	 * true if position components are within the allowed range, false otherwise.
     * 
     */

    return (position.pos[T_INDX] < lattice_param.n_T   && position.pos[T_INDX] >= 0 &&
            position.pos[X_INDX] < lattice_param.n_SPC && position.pos[X_INDX] >= 0 &&
            position.pos[Y_INDX] < lattice_param.n_SPC && position.pos[Y_INDX] >= 0 &&
            position.pos[Z_INDX] < lattice_param.n_SPC && position.pos[Z_INDX] >= 0   );
           
}

/* Checks if all position components and a Lorentz index are in the allowed range
   between 0 and N_mu for position and 0 and DIM for mu. */
inline bool positionmuValidQ(PosVec position, 
                             LorentzIdx mu) {

    /*
	 * Calls:
	 * =====
     * validPositionQ.
     * 
     * Macros:
	 * ======
     * X_INDX, Y_INDX, Z_INDX, T_INDX, DIM.
     * 
     * Global Variables:
     * ================
     * lattice_param
     *  
	 * Parameters:
	 * ==========
     * PosVec position:     position to be checked,
     * LorentzIdx mu:       Lorentz index to be checked.
     * 
	 * Returns:
	 * =======
	 * true if position components and Lorentz index are within the allowed range,
     * false otherwise.
     * 
     */

    return (validPositionQ(position) &&
            (mu == T_INDX || mu == X_INDX || 
             mu == Y_INDX || mu == Z_INDX   ));
}

/* Assigns numbers to a PosVec position struct, implementing the periodic boundary
   conditions. */
PosVec assignPosition(const PosIndex x, 
                      const PosIndex z, 
                      const PosIndex y, 
                      const PosIndex t) {

    /*
	 * Calls:
	 * =====
     * validPositionQ, makePeriodicBound.
     * 
     * Macros:
	 * ======
     * X_INDX, Y_INDX, Z_INDX, T_INDX.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
     * PosIndex x:      x component of the position,
     * PosIndex y:      y component of the position,
     * PosIndex z:      z component of the position,
     * PosIndex t:      t component of the position.
     * 
	 * Returns:
	 * =======
     * 
     */

    PosVec position;

    position.pos[X_INDX] = x;
    position.pos[Y_INDX] = y;
    position.pos[Z_INDX] = z;
    position.pos[T_INDX] = t;
    
    if(!validPositionQ(position)) {
        position = makePeriodicBound(position);
    }   /* if position not within bounds, put it within bounds */
    return position;
}

/* Prints a position to the screen. */
void printPosVec(const PosVec position) {

    /* 
	 * Calls:
	 * =====
     * printf.
     * 
     * Macros:
	 * ======
     * X_INDX, Y_INDX, Z_INDX, T_INDX, DIM.
     * 
     * Global Variables:
     * ================
     *  
	 * Parameters:
	 * ==========
     * PosVec position:     position to be printed,
     * 
	 * Returns:
	 * =======
     * 
     */

    printf("x: %hu y: %hu z: %hu t: %hu\n", position.pos[X_INDX],
                                            position.pos[Y_INDX],
                                            position.pos[Z_INDX],
                                            position.pos[T_INDX] );

}