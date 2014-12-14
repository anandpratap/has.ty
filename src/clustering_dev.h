	#ifndef _CLUSTERING_GPU_H
	#define _CLUSTERING_GPU_H
	
	__global__ void clustering_update(int ndata, int ncluster, double *x, unsigned int *np_dev, unsigned int *cidx_dev);
	
	#endif
