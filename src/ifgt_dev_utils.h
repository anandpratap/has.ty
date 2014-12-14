	#ifndef _IFGT_GPU_H
	#define _IFGT_GPU_H

	#include "common.h" // DEBUG 

	__device__ void AtomicAddDouble (volatile double *address, double value);
	__global__ void calc_source_monomials_gpu(int n, int k, double h, double *x, double *q);
	__global__ void eval_ifgt_gpu(int ndata, int ncluster, double h, double *y, double *f);

	#endif