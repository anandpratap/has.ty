	#ifndef _IFGT_GPU_H
	#define _IFGT_GPU_H

	#include <stdlib.h>

	#include <limits>  // max int
	#include <cmath> // exp
	
	// local includes	
	#include "utils_gpu.h"
	#include "common.h" // DEBUG 
	
	
	// calculate monomianls
	// used both for source and target
	__device__ void calc_monomial_gpu(int d, double *dx, int pmax, double *source_monomial);


	// calculate coefficients
	__device__ void calc_constant_gpu(int d, int pmax, double *constant);


	__device__ void calc_c_alpha_gpu(int n, int d, int pmax, double *x, double *c, double h, double *q, double *c_alpha);


	__global__ void calc_ifgt_gauss_transform_gpu(int n, int m, int d, int pmax, double *x, double *y, double *c, double h, double *q, double *f);
	
	#endif