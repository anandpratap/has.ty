	#ifndef _CLUSTERING_GPU_H
	#define _CLUSTERING_GPU_H
	#include <stdlib.h>
	#include <stdio.h>
	#include <cassert>
	#include "ifgt_dev.h"
	#include "utils_dev.h"
	#include "common.h"

	extern __device__ __constant__ double centers_dev[KMAX*DIMENSIONS];

	__device__ unsigned int get_minimum_idx_gpu(int nk, double *distance){
		unsigned int idx;
		double min_distance = MAX_DOUBLE;
		
		for(int i=0; i<nk; i++){
			if(distance[i] < min_distance){
				min_distance = distance[i];
				idx = i;
			}
		}
		assert(idx < nk);
		return idx;
	}

	
	__global__ void clustering_update(int ndata, int ncluster, double *x, unsigned int *np_dev, unsigned int *cidx_dev){

		int idx = blockDim.x*blockIdx.x + threadIdx.x;
		double dx[DIMENSIONS] = {0};
		double dx2[KMAX] = {MAX_DOUBLE};
		
		if(idx < ndata){
			for(int k=0; k<ncluster; k++){
				for(int d=0; d<DIMENSIONS; d++){
					dx[d] = x[idx*DIMENSIONS+d] - centers_dev[k*DIMENSIONS+d];
				}
				dx2[k] = l2normsq_dev(dx);
				assert(isfinite(dx2[k]));
			}
			cidx_dev[idx] = get_minimum_idx_gpu(ncluster, dx2);
		}
	}

	#endif
