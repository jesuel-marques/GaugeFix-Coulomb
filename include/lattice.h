#ifndef LATTICE_H
#define LATTICE_H
#include <stdbool.h>

#include <flags.h>

#include <types.h>


#define POSITION_IS_ODD(position)    (((position.i) ^ \
                                       (position.j) ^ \
                                       (position.k) ^ \
                                       (position.t)) & 1)

#define POSITION_IS_EVEN(position)   !(POSITION_IS_ODD(position))

#define LOOP_TEMPORAL(t)       for (t = 0; t < lattice_param.n_T; t++) 
#define LOOP_SPATIAL_DIR(mu)   for (position.mu = 0; \
                                    position.mu < lattice_param.n_SPC; \
                                    position.mu++) 
#define LOOP_SPATIAL(position) LOOP_SPATIAL_DIR(k)  \
                               LOOP_SPATIAL_DIR(j)  \
                               LOOP_SPATIAL_DIR(i)

// Paralelizing by slicing the time extent
#define LOOP_TEMPORAL_PARALLEL(t) OMP_PARALLEL_FOR \                                
                               LOOP_TEMPORAL(t)                           

#define LOOP_LORENTZ(mu) for (mu = 0; mu < DIM; mu++)  

#define LOOP_LORENTZ_SPATIAL(mu) for (mu = 0; mu < DIM - 1 ; mu++)


//  if position is odd, then the XOR of the first bit of each element
//  of position must be 1. Take AND with 1 select this first bit. Take the NOT of 
//  the odd code, because want 1 for even and 0 for odd.

//  Odd position means that the sum of the coordinates is odd and equivalente for even


//  Associations between the numeric indices and the lorentz directions

#define DIM 4                   //   Space-time lattice dimension

#define X_INDX 0
#define Y_INDX 1
#define Z_INDX 2
#define T_INDX 3


//  Geometric types definitions

typedef unsigned short PosIndex;             /*	Position index  */

typedef struct {
    PosIndex i, j, k;
    PosIndex t;
} PosVec;                                   /*	Struct for position vectors */


typedef unsigned short LorentzIdx;          //	Lorentz index

typedef enum {REAR, FRONT} Direction;       /*	Direction for link. A front link is the
                                                link in the positive direction from a 
                                                given point, whereas a rear link is a 
                                                link in the negative direction  */


typedef struct {

    short n_SPC;	//   Spatial lattice size
    short n_T;		//   Temporal lattice size

    int volume;		//	Number of sites in the lattice
    int spatial_volume;	

    int amount_of_links;
    int amount_of_points;

    bool error;
} GeometricParameters;

extern GeometricParameters lattice_param;

int initGeometricParameters(const short n_SPC, const short n_T);

bool validGeometricParametersQ(void);

inline bool validPositionQ(PosVec position);

bool positionmuValidQ(PosVec position, 
                       LorentzIdx mu);

PosVec makePositionValid(PosVec position);

PosVec assignPosition(const PosIndex x, 
                       const PosIndex y, 
                       const PosIndex z, 
                       const PosIndex t);

void printPosVec(const PosVec u);

PosVec hopPosPlus (const PosVec u, 
                     const LorentzIdx mu),
       hopPosMinus(const PosVec u, 
                     const LorentzIdx mu);

#endif  //LATTICE_H