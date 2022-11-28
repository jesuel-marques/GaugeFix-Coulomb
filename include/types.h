#ifndef TYPES_H
#define TYPES_H

#include <SU3_parameters.h>
#include <tgmath.h>

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
//  Data types definitions

typedef complex double WorkScalarType;      //  Calculations internally done in double
typedef WorkScalarType Scalar;              //  Default scalar is a scalar of type Work

// 2x2 types

typedef unsigned short MtrxIdx2;            //  2x2 matrix index

typedef struct {
    Scalar m[2 * 2];
} Mtrx2x2;                                  /*  Struct for 2x2 matrices */

typedef struct {
    double m[4];
} Mtrx2x2CK;                                /*  Struct for 2x2 matrices in the Cayley-
                                                Klein form 
                                                u = m[0] SU2_identity 
                                                    + i sum_i=1^3 m[i]sigma[i] */

// 3x3 types

typedef unsigned short MtrxIdx3;            //  3x3 matrix index
typedef unsigned short SU3AlgIdx;           //  SU(3) algebra index

typedef struct {
   complex double m[Nc * Nc];
} Mtrx3x3Dbl;                               /*  Struct for 3x3 matrices in double   */

typedef struct {
    complex float m[Nc * Nc];
} Mtrx3x3Flt;                               /*  Struct for 3x3 matrices in float    */

typedef struct {
   Scalar m[Nc];
} Vec3;                                     /*  Struct for 3 component vectors in color
                                                space. */

typedef struct {
    Scalar m[(Nc * Nc - 1) + 1];
} MtrxSU3Alg;                               /*  Struct for SU(3) algebra components. */

typedef Mtrx3x3Dbl WorkMtrx;                //  Calculations internally done in double
typedef WorkMtrx Mtrx3x3;                   //  Default matrix is matrix of type Work

typedef enum {R, S, T} Submtrx;             /*  Submtrx provides an alias for the 
                                                indices of Cabbibo-Marinari subdivision 
                                                of an SU(3) matrix.

/*
This is a way of embedding SU(2) matrices in an SU(3) matrix. For each of the indices,
the X marks where the SU(2) matrix goes and the rest is kept as either a 0 or a 1, 
as shown.

R ->  X X 0
      X X 0
      0 0 1

S ->  X 0 X
      0 1 0
      X 0 X

T ->  1 0 0
      0 X X
      0 X X

This is used in this program for the updates to the local gauge transformations. The 
updates are performed for each submatrix in turn.
*/


#endif  //TYPES_H