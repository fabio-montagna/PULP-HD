#ifndef AUX_FUNCTIONS_H_
#define AUX_FUNCTIONS_H_
 
#include "init.h"
#include <math.h>


#define ROUND(num) ( (num - floorf(num) > 0.5f) ? ceilf(num) : floorf(num) )


void hamming_dist(uint32_t q[bit_dim + 1], uint32_t aM[][bit_dim + 1], int sims[classes]);
int max_dist_hamm(int distances[classes]);
void computeNgram(float buffer[channels], uint32_t iM[][bit_dim + 1], uint32_t ciM[][bit_dim + 1], uint32_t query[bit_dim + 1]);
int numberOfSetBits(uint32_t i);

#endif
