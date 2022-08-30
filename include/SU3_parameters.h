//  Simulation Parameters

#define N_SPC   24  //   Spatial lattice size
#define N_T     16    // Temporal lattice size
#define DIM     4     // Space-time lattice dimension

#define SPATIAL_VOLUME  N_SPC * N_SPC * N_SPC
#define VOLUME          N_SPC * N_SPC * N_SPC * N_T  //	Number of points in the lattice

#define AMOUNT_OF_LINKS  DIM*VOLUME
#define AMOUNT_OF_POINTS VOLUME

#define Nc 3    //  Number of colors


// Other parameters

#define FIRST_CONFIG    500
#define LAST_CONFIG     11000
#define CONFIG_STEP     10

#define MAX_CONFIGS (LAST_CONFIG - FIRST_CONFIG) / CONFIG_STEP

#define MAX_LENGTH_NAME 2000

#define NUM_THREADS 8


#define     NEED_BYTE_SWAP_IN

// #define  NEED_BYTE_SWAP_OUT


#define     CONV_CFG_TO_WORKING_PRECISION

// #define     CONV_CFG_FROM_WORKING_PRECISION

// #define  CONV_GT_TO_WORKING_PRECISION

// #define  CONV_GT_FROM_WORKING_PRECISION


// #define CHECK_POSITION_BOUNDS