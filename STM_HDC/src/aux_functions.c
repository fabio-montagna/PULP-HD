
#include <stdlib.h>
#include "aux_functions.h"
#include <math.h>


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
/*************************************************************************
	DESCRIPTION: computes the Hamming Distance for each class.

	INPUTS:
		q        : query hypervector
		aM		 : Associative Memory matrix

	OUTPUTS:
		sims	 : Distances' vector
**************************************************************************/

	int r_tmp = 0;
	uint32_t tmp = 0;

	for(int i = 0; i < classes; i++){

		for(int j = 0; j < bit_dim + 1; j++){

			tmp = q[j] ^ aM[i][j];


			r_tmp += numberOfSetBits(tmp);

		}

		sims[i] = r_tmp;
		r_tmp = 0;

	}


}


void computeNgram(float buffer[channels], uint32_t iM[][bit_dim + 1], uint32_t ciM[][bit_dim + 1], uint32_t query[bit_dim + 1]){
/*************************************************************************
	DESCRIPTION: computes the N-gram

	INPUTS:
		buffer   :  data input
		iM       :	Item Memory for the IDs of channels
		chAM     :  Continuous Item Memory for the values of a channel
	OUTPUTS:
		query    :  query hypervector
**************************************************************************/

	int r[channels];
	int ix;
	uint32_t tmp = 0;
	int i, j;
	uint32_t chHV[channels + 1][bit_dim + 1] = {0};


	//Quantization: each sample is rounded to the nearest integer
	for(i = 0; i < channels ; i++){

		r[i] = (int)(ROUND((float)buffer[i]));

	}

	//Spatial Encoder: captures the spatial information for a given time-aligned samples of channels
	for(i = 0; i < bit_dim + 1; i++){

		query[i] = 0;
		for(j = 0; j < channels; j++){

			ix = r[j];
			tmp = ciM[ix][i] ^ iM[j][i];
			chHV[j][i] = tmp;

		}
		//this is done to make the dimension of the matrix odd for the componentwise majority.
		chHV[channels][i] = chHV[0][i] ^ chHV[1][i];



		uint32_t majority = 0;
		//componentwise majority: inserts the value of the i-th bit of each chHV row in the variable "majority"
		//and then computes the number of 1's with the function numberOfSetBits(uint32_t).
		for(int z = 31; z >= 0; z--){

			for(int j = 0 ; j < channels + 1; j++){

				majority = majority | (((chHV[j][i] & ( 1 << z)) >> z) << j);


			}
			if (numberOfSetBits(majority) > 2) query[i] = query[i] | ( 1 << z ) ;

			majority = 0;
		}



	}

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
