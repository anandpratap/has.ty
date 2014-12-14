	#ifndef _UTILS_GPU_H
	#define _UTILS_GPU_H
	#include <stdlib.h>

	#include <cmath> // exp
	#include <stdlib.h>
	#include <stdio.h>
	// local includes
	#include "common.h"
		
	
	__device__ int nchoosek_dev(int n, int k){
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

	
	__device__ double l2normsq_dev(double *x){
		double norm = 0.0;
		for(int i=0; i<DIMENSIONS; i++){
			norm += x[i]*x[i];
		}
		return norm;
	}
	

	
	void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
	{
		if (code != cudaSuccess) 
		{
			fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
			if (abort) exit(code);
		}
	}

	#endif