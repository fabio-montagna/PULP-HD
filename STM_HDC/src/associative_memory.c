#include "associative_memory.h"
#include "aux_functions.h"


int associative_memory_32bit(uint32_t q_32[bit_dim + 1], uint32_t aM_32[][bit_dim + 1]){

/*************************************************************************
	DESCRIPTION:  tests the accuracy based on an input testing queries

	INPUTS:
		q_32        : query hypervector
		aM_32		: Trained Associative Memory
	OUYTPUTS:
		class       : classification result
**************************************************************************/

	int sims[classes] = {0};
	int class;


	//HAMMING DISTANCES
	hamming_dist(q_32, aM_32, sims);

	//CLASSIFICATION WITH HAMMING METRIC
	class = max_dist_hamm(sims);

    return class;
}


