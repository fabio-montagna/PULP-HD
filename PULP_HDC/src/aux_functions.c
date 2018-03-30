#include "aux_functions.h"
#include "rt/rt_alloc.h" 

#define ROUND(num) ((num - floorf(num) > 0.5f) ? ceilf(num) : floorf(num))


int max_dist_hamm(int distances[classes]){
/*************************************************************************
	DESCRIPTION: computes the maximum Hamming Distance.

	INPUTS:
		distances     : distances associated to each class
	OUTPUTS:
		max_index     : the class related to the maximum distance
**************************************************************************/
	int max = distances[0];
	int max_index = 0;

	for(int i = 1; i < classes; i++){

		if(max > distances[i]){

			max = distances[i];
			max_index = i;

		}

	}

	return max_index;
}



void hamming_dist(uint32_t q[bit_dim + 1], uint32_t aM[][bit_dim + 1], int sims[classes]){
/**************************************************************************
	DESCRIPTION: computes the Hamming Distance for each class.

	INPUTS:
		q        : query hypervector
		aM		 : Associative Memory matrix

	OUTPUTS:
		sims	 : Distances' vector
***************************************************************************/

	int r_tmp = 0;

	#pragma omp parallel num_threads(CORE)
	{
	uint32_t tmp = 0;
	for(int i = 0; i < classes; i++){
	#pragma omp for reduction(+:r_tmp)
		for(int j = 0; j < bit_dim + 1; j++){

			tmp = q[j] ^ aM[i][j];

#if WOLF
			r_tmp +=__builtin_popcount(tmp);
#else
			r_tmp += numberOfSetBits(tmp);
#endif
		}
	#pragma omp master
	{
		sims[i] = r_tmp;
		r_tmp = 0;
	}
	#pragma omp barrier
	}
	}//omp

}

PULP_L1_DATA uint32_t chHV[channels + 1][bit_dim + 1] = {0};

void computeNgram(float buffer[channels], uint32_t iM[][bit_dim + 1], uint32_t chAM[][bit_dim + 1], uint32_t query[bit_dim + 1]){
/*************************************************************************
	DESCRIPTION: computes the N-gram

	INPUTS:
		buffer   :  input data
		iM       :	Item Memory for the IDs of channels
		chAM     :  Continuous Item Memory for the values of a channel
	OUTPUTS:
		query    :  query hypervector
**************************************************************************/

	int r[channels];


	#pragma omp parallel num_threads(CORE)
	{

	int ix;
	uint32_t tmp = 0;
	int i, j;
	uint32_t chHV[channels + 1][bit_dim + 1] = {0};

	#pragma omp master
	{
	//Quantization: each sample is rounded to the nearest integer
	for(i = 0; i < channels ; i++){

		r[i] = (int)(ROUND((float)buffer[i]));

	}



	}//master
	#pragma omp barrier

	//Spatial Encoder: captures the spatial information for a given time-aligned samples of channels
	#pragma omp for
	for(i = 0; i < bit_dim + 1; i++){

		query[i] = 0;
		for(j = 0; j < channels; j++){

			ix = r[j];
			tmp = iM[ix][i] ^ chAM[j][i];
			chHV[j][i] = tmp;

		}
		//this is done to make the dimension of the matrix for the componentwise majority odd.
		chHV[channels][i] = chHV[0][i] ^ chHV[1][i];

		//componentwise majority: insert the value of the ith bit of each chHV row in the variable "majority"
		//and then compute the number of 1's with the function numberOfSetBits(uint32_t).
#if WOLF == 1

	uint32_t majority = 0;
	uint32_t res0, res1, res2, res3, res4;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 31);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 31);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 31);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 31);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 31);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 31);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 30);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 30);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 30);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 30);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 30);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 30);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 29);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 29);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 29);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 29);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 29);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 29);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 28);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 28);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 28);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 28);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 28);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 28);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 27);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 27);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 27);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 27);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 27);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 27);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 26);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 26);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 26);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 26);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 26);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 26);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 25);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 25);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 25);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 25);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 25);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 25);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 24);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 24);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 24);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 24);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 24);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 24);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 23);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 23);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 23);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 23);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 23);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 23);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 22);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 22);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 22);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 22);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 22);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 22);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 21);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 21);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 21);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 21);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 21);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 21);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 20);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 20);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 20);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 20);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 20);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 20);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 19);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 19);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 19);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 19);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 19);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 19);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 18);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 18);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 18);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 18);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 18);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 18);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 17);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 17);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 17);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 17);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 17);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 17);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 16);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 16);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 16);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 16);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 16);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 16);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 15);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 15);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 15);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 15);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 15);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 15);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 14);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 14);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 14);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 14);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 14);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 14);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 13);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 13);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 13);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 13);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 13);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 13);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 12);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 12);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 12);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 12);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 12);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 12);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 11);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 11);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 11);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 11);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 11);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 11);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 10);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 10);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 10);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 10);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 10);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 10);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 9);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 9);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 9);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 9);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 9);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 9);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 8);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 8);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 8);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 8);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 8);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 8);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 7);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 7);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 7);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 7);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 7);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 7);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 6);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 6);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 6);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 6);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 6);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 6);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 5);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 5);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 5);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 5);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 5);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 5);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 4);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 4);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 4);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 4);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 4);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 4);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 3);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 3);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 3);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 3);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 3);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 3);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 2);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 2);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 2);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 2);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 2);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 2);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 1);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 1);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 1);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 1);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 1);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 1);
	majority = 0;

	res0 = __builtin_bitextractu(chHV[0][i], 1, 0);
	res1 = __builtin_bitextractu(chHV[1][i], 1, 0);
	res2 = __builtin_bitextractu(chHV[2][i], 1, 0);
	res3 = __builtin_bitextractu(chHV[3][i], 1, 0);
	res4 = __builtin_bitextractu(chHV[4][i], 1, 0);
	majority = __builtin_bitinsert(majority, res0, 1, 0);
	majority = __builtin_bitinsert(majority, res1, 1, 1);
	majority = __builtin_bitinsert(majority, res2, 1, 2);
	majority = __builtin_bitinsert(majority, res3, 1, 3);
	majority = __builtin_bitinsert(majority, res4, 1, 4);

	if (__builtin_popcount(majority) > 2) query[i] = __builtin_bitinsert(query[i], 1, 1, 0);
	majority = 0;



#else
		uint32_t majority = 0;

		for(int z = 31; z >= 0; z--){

			for(int j = 0 ; j < channels + 1; j++){

				majority = majority | (((chHV[j][i] & ( 1 << z)) >> z) << j);


			}
			if (numberOfSetBits(majority) > 2) query[i] = query[i] | ( 1 << z ) ;

			majority = 0;
		}

#endif
	}

}//omp

}

int numberOfSetBits(uint32_t i)
{
/*************************************************************************
	DESCRIPTION:   computes the number of 1's

	INPUTS:
		i        :  the i-th variable that composes the hypervector

**************************************************************************/

     i = i - ((i >> 1) & 0x55555555);
     i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
     return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}
