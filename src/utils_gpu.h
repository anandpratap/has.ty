	#ifndef _UTILS_GPU_H
	#define _UTILS_GPU_H
	#include <stdlib.h>

	#include <cmath> // exp
	// local includes
	#include "common.h"
		
	

	__device__ int nchoosek_gpu(int n, int k);
	
	__device__ double l2normsq_gpu(double *x);
	
	void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true);
	#define GPU_ERROR_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
	#endif