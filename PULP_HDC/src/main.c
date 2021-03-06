#include <stdio.h>
#include <stdint.h>
#include "associative_memory.h"
#include "aux_functions.h"
#include "init.h"
//the data.h file can be created directly in MATLAB (after the simulation)
//using the function "data_file_creator.m"
#include "data.h"

//THE NUMBER OF CORES (1, 2, 4 OR 8) USED FOR THE EXECUTION CAN BE SET IN THE Makefile : CORE ?= 8.   

PULP_L1_DATA float buffer[4];

int main(){
 
          
	uint32_t overflow = 0;
	uint32_t old_overflow = 0;
	uint32_t mask = 1;
	uint32_t q[bit_dim + 1] = {0};
	uint32_t q_N[bit_dim + 1] = {0};
	int class;
  
	for(int ix = 0; ix < NUMBER_OF_INPUT_SAMPLES; ix = ix + N){
#if HWPERF
	rt_perf_t perf;
	rt_perf_init(&perf);
	rt_perf_reset(&perf);
 
   	rt_perf_start(&perf);
#endif
		for(int z = 0; z < N; z++){
 
			for(int j = 0; j < channels; j++){
			 
				buffer[j] = TEST_SET[j][ix + z];
			 
			}
			//Spatial and Temporal Encoder: computes the N-gram.
			//N.B. if N = 1 we don't have the Temporal Encoder but only the Spatial Encoder.
			if (z == 0) computeNgram(buffer, iM, chAM, q);
			else{

 
				computeNgram(buffer, iM, chAM, q_N);
 
				//Here the hypervetor q is shifted by 1 position as permutation,
				//before performing the componentwise XOR operation with the new query (q_N).
				overflow = q[0] & mask;

				for(int i = 1; i < bit_dim; i++){

					old_overflow = overflow;
					overflow = q[i] & mask;
					q[i] = (q[i] >> 1) | (old_overflow << (32 - 1));
					q[i] = q_N[i] ^ q[i];

				}

				old_overflow = overflow;
				overflow = (q[bit_dim] >> 16) & mask;
				q[bit_dim] = (q[bit_dim] >> 1) | (old_overflow << (32 - 1));
				q[bit_dim] = q_N[bit_dim] ^ q[bit_dim];

				q[0] = (q[0] >> 1) | (overflow << (32 - 1));
				q[0] = q_N[0] ^ q[0];


			}
 
		}
	
 
	    //classifies the new N-gram through the Associative Memory matrix.
		if (N == 1)
			class = associative_memory_32bit(q, aM_32);
		else
			class = associative_memory_32bit(q, aM_32);
#if HWPERF
	rt_perf_stop(&perf);
    	rt_perf_save(&perf);
   	printf("%d\t",rt_perf_read (RT_PERF_CYCLES));
#endif
	
 	printf("%d\n", class);



	}

    return 0; 
}

