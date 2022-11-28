#ifndef TYPES_H
#define TYPES_H

#include <SU3_parameters.h>
#include <tgmath.h>

//  Geometric types definitions

typedef unsigned short PosIndex;

typedef struct {
    PosIndex i, j, k;
    PosIndex t;
} PosVec;  //	struct for position vectors


typedef unsigned short LorentzIdx;

typedef enum {REAR, FRONT} Direction;

//  Data types definitions

typedef complex double WorkScalarType;
typedef WorkScalarType Scalar;

// 2x2 types

typedef unsigned short SU2ColorIdx;

typedef struct {
    Scalar m[2 * 2];
} Mtrx2x2;

typedef struct {
    double m[4];
} Mtrx2x2CK;

// 3x3 types

typedef unsigned short MtrxIdx3;
typedef unsigned short SU3AlgIdx;

typedef struct {
   complex double m[Nc * Nc];
} Mtrx3x3Dbl; 

typedef struct {
    complex float m[Nc * Nc];
} Mtrx3x3Flt;

typedef struct {
   Scalar m[Nc];
} Vec3;

typedef struct {
    Scalar m[(Nc * Nc - 1) + 1];
} MtrxSU3Alg;

typedef Mtrx3x3Dbl WorkMtrx;
typedef WorkMtrx Mtrx3x3;

typedef enum {R, S, T} Submtrx;

#endif  //TYPES_H