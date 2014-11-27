	#ifndef _UTILS_GPU_H
	#define _UTILS_GPU_H
	#include <stdlib.h>

	#include <cmath> // exp
	// local includes
	#include "common.h"
		
	__device__ void set_zeros_gpu(int n, double *x){
		for(int i=0; i< n; i++)
			x[i] = 0.0;
	}


	__device__ int nchoosek_gpu(int n, int k){
		int n_k = n - k;
		if (k < n_k){
			k = n_k;
			n_k = n - k;
		}
		int  nchsk = 1; 
		for ( int i = 1; i <= n_k; i++){
			nchsk *= (++k);
			nchsk /= i;
		}

		return nchsk;
	}

	
	__device__ double l2normsq_gpu(int n, double *x){
		double norm = 0.0;
		for(int i=0; i<n; i++){
			norm += x[i]*x[i];
		}
		return norm;
	}

	#endif