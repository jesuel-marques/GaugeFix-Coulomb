#include <stdio.h>
#include <stdbool.h>

#include <flags.h>
#include <lattice.h>
#include <SU3_ops.h>
#include <types.h>

GeometricParameters lattice_param;

PosVec *neighbour_table;


static inline int mod(short a, unsigned short b) {
    
    /*
	 * Description:
     * ===========
	 * Modulo function which deals with negative numbers.
     *  
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


PosVec makePeriodicBound(PosVec position) {

     /*
	 * Description:
     * ===========
	 * Implements the periodic boundary conditions at the geometric level.
     *  
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

static void makeNeighbourTable() {

    /*
	 * Description:
     * ===========
	 * Sets up a table which records which lattice points are neighbours to which other.
     *  
	 * Calls:
	 * =====
     * malloc,
     * makePeriodicBound.
     * 
     * Macros:
	 * ======
     * X_INDX, Y_INDX, Z_INDX, T_INDX, 
     * LOOP_TEMPORAL_PARALLEL, LOOP_SPATIAL, LOOP_LORENTZ, DIM.
     * 
     * Global Variables:
     * ================
     * lattice_param, neighbour_table.
     *  
	 * Parameters:
	 * ==========
     * 
	 * Returns:
	 * =======
	 * 
     */

    neighbour_table = (PosVec *) malloc(lattice_param.amount_of_points
                                        * DIM * 2 * sizeof(PosVec));
    short t;
    LOOP_TEMPORAL_PARALLEL(t) {
        PosVec position = {.pos={0, 0, 0, t}};
        PosVec neighbour = {.pos={0, 0, 0, 0}};

        LorentzIdx mu;
        LOOP_SPATIAL(position) {
            LOOP_LORENTZ(mu) {
                neighbour = position;
                neighbour.pos[mu]--;

                *(neighbour_table + ((((position.pos[T_INDX]  * lattice_param.n_SPC 
                                      + position.pos[Z_INDX]) * lattice_param.n_SPC
                                      + position.pos[Y_INDX]) * lattice_param.n_SPC
                                      + position.pos[X_INDX]) * DIM 
                                      + mu) * 2 + REAR) = makePeriodicBound(neighbour);
            
                neighbour = position;
                neighbour.pos[mu]++;

                *(neighbour_table + ((((position.pos[T_INDX]  * lattice_param.n_SPC 
                                      + position.pos[Z_INDX]) * lattice_param.n_SPC
                                      + position.pos[Y_INDX]) * lattice_param.n_SPC
                                      + position.pos[X_INDX]) * DIM 
                                      + mu) * 2 + FRONT) = makePeriodicBound(neighbour);
            }
        }
    }

}


inline PosVec getNeighbour(PosVec position, LorentzIdx mu, Direction dir) {

    /*
	 * Description:
     * ===========
	 * Returns the front or rear (depending on dir) neighbour of position in Lorentz 
     * direction mu. 
     * 
	 * Calls:
	 * =====
     * printf.
     * 
     * Macros:
	 * ======
     * DIM.
     * 
     * Global Variables:
     * ================
     * neighbour_table
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

    /* PosVec neighbour = position; 

    if(dir == REAR) {
        neighbour.pos[mu] --;
    }
    else{
        neighbour.pos[mu]++;
    } 

    return makePeriodicBound(neighbour); */
    
    return *(neighbour_table + ((((  position.pos[T_INDX]  * lattice_param.n_SPC 
                                   + position.pos[Z_INDX]) * lattice_param.n_SPC
                                   + position.pos[Y_INDX]) * lattice_param.n_SPC
                                   + position.pos[X_INDX]) * DIM + mu) * 2 + dir);
}


void finalizeGeometry(void) {

    /*
	 * Description:
     * ===========
	 * Frees up the memory allocated for the neighbour table.
     *  
	 * Calls:
	 * =====
     * free.
     * 
     * Macros:
	 * ======
     * 
     * Global Variables:
     * ================
     * neighbour_table
     *  
	 * Parameters:
	 * ==========
     * 
	 * Returns:
	 * =======
	 * 
     */

    free(neighbour_table);
}


bool validGeometricParametersQ(void) {

    /*
	 * Description:
     * ===========
	 * Check if current geometric parameters are valid.
     *  
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

    if(lattice_param.n_SPC <= 0                                 || 
       lattice_param.n_T   <= 0                                 ||
       lattice_param.error == 1                                 ||
       lattice_param.spatial_volume != ns * ns * ns             ||
       lattice_param.volume != ns * ns * ns * nt                ||
       lattice_param.amount_of_links != DIM * ns * ns * ns * nt ||
       lattice_param.amount_of_points != ns * ns * ns * nt        ) {
        
        lattice_param.error = 1;
        return false;
    }
    else{
        return true;
    }	
}


int initGeometry(const short n_s, const short n_t) {

    /*
	 * Description:
     * ===========
	 * Initializes the global variable lattice_param which contains the lattice 
     * geometric parameters, performing check to see if the parameters set are
     * reasonable. After that it generates a table of neighbours for neighbour lookups.
     *  
	 * Calls:
	 * =====
     * validGeometricParametersQ, makeNeighbourTable.
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

	lattice_param.amount_of_links  = DIM * lattice_param.volume;
	lattice_param.amount_of_points = lattice_param.volume;

    if(!validGeometricParametersQ()) {
        return 1;
    }

    makeNeighbourTable();

    return 0;
}


inline bool validPositionQ(PosVec position) {

    /*
	 * Description:
     * ===========
	 * Checks if all position components are in the allowed range between 0 and N_mu.
     *  
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

    LorentzIdx mu;
    LOOP_LORENTZ(mu) {
        if(position.pos[mu] >= (mu != T_INDX ? 
                               lattice_param.n_SPC : 
                               lattice_param.n_T) &&
           position.pos[mu] < 0 )
            return false;
    }

    return true;
}


inline bool positionmuValidQ(PosVec position, 
                             LorentzIdx mu) {

    /*
	 * Description:
     * ===========
	 * Checks if all position components and a Lorentz index are in the allowed range
     * between 0 and N_mu for position and 0 and DIM for mu.
     * 
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

    if(validPositionQ(position) &&
        (mu == T_INDX || mu == X_INDX || 
         mu == Y_INDX || mu == Z_INDX   )) {
        return true;
    }

    return false;
}


PosVec assignPosition(const PosIndex x, 
                      const PosIndex z, 
                      const PosIndex y, 
                      const PosIndex t) {

    /*
	 * Description:
     * ===========
	 * Assigns numbers to a PosVec position struct, implementing the periodic boundary
     * conditions.
     * 
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

void printPosVec(const PosVec position) {

    /*
	 * Description:
     * ===========
	 * Prints a position to the screen.
     * 
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