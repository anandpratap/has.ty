	#ifndef _UTILS_GPU_H
	#define _UTILS_GPU_H
	#include <stdlib.h>

	#include <cmath> // exp
	// local includes
	#include "common.h"
		
	__device__ void set_zeros_gpu(int n, double *x);


	__device__ int nchoosek_gpu(int n, int k);
	
	__device__ double l2normsq_gpu(int n, double *x);
	
	#endif