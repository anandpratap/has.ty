	#ifndef _UTILS_DEV_H
	#define _UTILS_DEV_H
	#include <stdlib.h>

	#include <cmath> // exp
	// local includes
	#include "common.h"
		
	
	/**
	* @brief Return nchoosek  on the device
	* 
	* @param[in] n
	* @param[in] k
	**/
	__device__ int nchoosek_dev(int n, int k);
	
	/**
	* @brief Return square of the L2 norm
	*
	* @param[in] *x 
	**/
	__device__ double l2normsq_dev(double *x);
	
	void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true);
	#define CUDA_CALL(ans) { gpuAssert((ans), __FILE__, __LINE__); }
	
	#endif