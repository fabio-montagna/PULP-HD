#ifndef INIT_H_
#define INIT_H_

#include <stdint.h>
#include <math.h>
#include "omp.h"

 

#define PULP_L1_DATA RT_L1_DATA
#define PULP_L2_DATA RT_L2_DATA
 

//dimension of the hypervectors
#define dimension 10000
//number of classes to be classify
#define classes 5
//number of acquisition's channels
#define channels 4
//dimension of the hypervectors after compression (dimension/32 rounded to the smallest integer)
#define bit_dim 312
//number of input samples
#define NUMBER_OF_INPUT_SAMPLES 14883
//dimension of the N-grams (models for N = 1 and N = 5 are contained in data.h)
#define N 5
//CHANNELS_VOTING for the componentwise majority must be odd
#define CHANNELS_VOTING channels + 1
//Set to 1 to use the profiling's functions
#define HWPERF 0
//Set to 1 to use the built-ins with Wolf, 0 for the execution without built-ins. 
#define WOLF 1

#define TEST 1



#endif
