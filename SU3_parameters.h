//  Simulation Parameters

#define N_SPC 24  //   Spatial lattice size
#define N_T 28    // Temporal lattice size
#define DIM 4     // Space-time lattice dimension

#define VOLUME N_SPC * N_SPC * N_SPC * N_T  //	Number of points in the lattice


#define Nc 3    //  Number of colors


// Other parameters

#define FIRST_CONFIG 1000

#define LAST_CONFIG 11000

#define CONFIG_STEP 10

#define MAX_CONFIGS (LAST_CONFIG - FIRST_CONFIG) / CONFIG_STEP

#define MAX_LENGTH_NAME 400

#define NUM_THREADS 4

#define NEED_BYTE_SWAP_IN

// #define NEED_BYTE_SWAP_OUT

#define NEED_CONV_TO_WORKING_PRECISION

#define NEED_CONV_FROM_WORKING_PRECISION