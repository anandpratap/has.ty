	#ifndef _IFGT_GPU_H
	#define _IFGT_GPU_H

	#include <stdlib.h>

	#include <limits>  // max int
	#include <cmath> // exp
	
	// local includes	
	#include "utils_gpu.h"
	#include "common.h" // DEBUG 

	extern __device__ __constant__ double constant_dev[NALPHA];
	extern __device__ __constant__ double centers_dev[KMAX*DIMENSIONS];

	__global__ void calc_source_monomials_gpu(int n, int k, double h, double *x, double *q);
	__global__ void eval_ifgt_gpu(int ndata, int ncluster, double h, double *y, double *f);

	#endif