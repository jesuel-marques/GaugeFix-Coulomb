//  Simulation Parameters

#define N_SPC 32  //   Spatial lattice size
#define N_T 16    // Temporal lattice size
#define DIM 4     // Space-time lattice dimension

#define VOLUME N_SPC * N_SPC * N_SPC * N_T  //	Number of points in the lattice


#define Nc 3    //  Number of colors


// Other parameters

#define FIRST_CONFIG 1000

#define LAST_CONFIG 6000

#define CONFIG_STEP 5

#define MAX_CONFIGS (LAST_CONFIG - FIRST_CONFIG) / CONFIG_STEP

#define MAX_LENGTH_NAME 2000

#define NUM_THREADS 4

#define NEED_BYTE_SWAP_IN

// #define NEED_BYTE_SWAP_OUT

#define NEED_CONV_TO_WORKING_PRECISION

#define NEED_CONV_FROM_WORKING_PRECISION