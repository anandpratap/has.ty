	#ifndef _IFGT_DEV_H
	#define _IFGT_DEV_H
	#include <iostream>
	#include <fstream>
	#include <stdio.h>
	#include <chrono>

	#include "cuda.h"
	#include "cuda_runtime_api.h"
	#include "curand.h"

	#include "common.h"
	#include "utils_host.h"
	#include "utils_dev.h"


	#include "clustering_host.h" // host clustering is current used
	#include "ifgt_host.h"
	#include "clustering_dev.h" // host clustering is current used
	#include "ifgt_dev_utils.h"
	
	// device constant memory variable
	__device__ __constant__ double constant_dev[NALPHA];
	__device__ __constant__ double centers_dev[KMAX*DIMENSIONS];

	void clustering_dev(int ndata, int ncluster, double *x, double *q, unsigned int *cidx, double *centers, unsigned int *np, struct cluster *clusters, int nthreads){
		
		double *x_dev;
		unsigned int *cidx_dev;
		unsigned int *np_dev;

		CUDA_CALL(cudaMalloc((void**) &x_dev, ndata*DIMENSIONS*sizeof(double)));
		CUDA_CALL(cudaMalloc((void**) &cidx_dev, ndata*sizeof(unsigned int)));
		CUDA_CALL(cudaMalloc((void**) &np_dev, ncluster*sizeof(unsigned int)));

		CUDA_CALL(cudaMemcpy(x_dev, x, ndata*DIMENSIONS*sizeof(double), cudaMemcpyHostToDevice));

		int nblocks = ceil(ndata/(float)nthreads); 
		dim3 threads(nthreads);
		dim3 blocks(nblocks);

		// initialize
		for(int i=0; i<ndata; i++){
			cidx[i] = rand() % ncluster;
		}

		calc_centers(ndata, ncluster, x, cidx, centers, np);
		// cuda clustering
		CUDA_CALL(cudaMemcpy(cidx_dev, cidx, ndata*sizeof(unsigned int), cudaMemcpyHostToDevice));
		CUDA_CALL(cudaMemcpy(np_dev, np, ncluster*sizeof(unsigned int), cudaMemcpyHostToDevice));
		cudaDeviceSynchronize();
		unsigned int iter = 0;
		unsigned int delta = 100000;
		while (1){
			CUDA_CALL(cudaMemcpyToSymbol(centers_dev, centers, ncluster*DIMENSIONS*sizeof(double)));
			cudaDeviceSynchronize();
			clustering_update<<<blocks, threads>>>(ndata, ncluster, x_dev, np_dev, cidx_dev);
			CUDA_CALL(cudaMemcpy(cidx, cidx_dev, ndata*sizeof(int), cudaMemcpyDeviceToHost));
			cudaDeviceSynchronize();
			calc_centers(ndata, ncluster, x, cidx, centers, np);
			iter += 1;
			cudaDeviceSynchronize();
			if(delta/(double)ndata < 1e-4 || iter > 300){
				printf("GPU Clustering Iterations: %d\n", iter);
				break;
			}

		}

		CUDA_CALL(cudaMemcpy(cidx, cidx_dev, ndata*sizeof(unsigned int), cudaMemcpyDeviceToHost));
		
		std::fstream outfile("clustering.txt", std::ios_base::out);
		for(int i=0; i<ndata; i++){
			outfile << x[i*DIMENSIONS] << " " << x[i*DIMENSIONS+1] << " " <<0<<" "<<0<<" "<< cidx[i] << "\n";
		}


		cudaFree(x_dev);
		cudaFree(cidx_dev);
		cudaFree(np_dev);
		
	}

	void ifgt_dev(int ndata, int mdata, int ncluster, double h, double *x, double *q, double *y, double *f, int nthreads){

		// cuda stuff
		hline();
		gprintf("GPU: Improved Fast Gauss Transform\n");
		hline();


		// temporarily using cpu clusters
		double *centers = new double[DIMENSIONS*ncluster]();
		unsigned int *cidx = new unsigned int [ndata]();
		unsigned int *np = new unsigned int [ncluster]();
		struct cluster *clusters = new struct cluster[ncluster]();
		clustering(ndata, ncluster, x, cidx, centers, np);
		rearrange_data(ndata, ncluster, x, q, cidx, centers, np, clusters);

		// Device clustering
		auto start = std::chrono::system_clock::now();
		
		unsigned int *tmp_cidx = new unsigned int [ndata]();
		double *tmp_centers = new double[DIMENSIONS*ncluster]();
		unsigned int *tmp_np = new unsigned int [ncluster]();
		struct cluster *tmp_clusters = new struct cluster[ncluster]();

		clustering_dev(ndata, ncluster, x, q, tmp_cidx, tmp_centers, tmp_np, tmp_clusters, nthreads);
		auto end = std::chrono::system_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << "GPU_clustering time (in ms): " << elapsed.count() << '\n';

		delete[] tmp_cidx;
		delete[] tmp_centers;
		delete[] tmp_np;
		delete[] tmp_clusters;

		unsigned int max_num_points = get_max_num_points(ncluster, clusters);

		double *constant_cpu = new double[NALPHA]();
		calc_constant(constant_cpu);

		double *x_dev;
		double *y_dev;
		double *q_dev;
		double *f_dev;


		printf("Max N point %d\n", max_num_points);
		CUDA_CALL(cudaMalloc((void**) &x_dev, max_num_points*DIMENSIONS*sizeof(double)));
		CUDA_CALL(cudaMalloc((void**) &centers_dev, ncluster*DIMENSIONS*sizeof(double)));
		CUDA_CALL(cudaMalloc((void**) &q_dev, max_num_points*sizeof(double)));
		CUDA_CALL(cudaMalloc((void**) &f_dev, mdata*sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(constant_dev, constant_cpu, NALPHA*sizeof(double)));
		CUDA_CALL(cudaMemcpyToSymbol(centers_dev, centers, ncluster*DIMENSIONS*sizeof(double)));

			
		start = std::chrono::system_clock::now();

		for(int c=0; c<ncluster; c++){
			struct cluster current_cluster = clusters[c];
			unsigned int n = current_cluster.n;

			CUDA_CALL(cudaMemcpy(x_dev, current_cluster.x, n*DIMENSIONS*sizeof(double), cudaMemcpyHostToDevice));
			CUDA_CALL(cudaMemcpy(q_dev, current_cluster.q, n*sizeof(double), cudaMemcpyHostToDevice));
			cudaDeviceSynchronize();
			if(current_cluster.n != 0){
				int nblocks_s = ceil((current_cluster.n)/(float)nthreads); 
				dim3 threads_s(nthreads);
				dim3 blocks_s(nblocks_s);
				calc_source_monomials_gpu<<<blocks_s,threads_s, nthreads*sizeof(double)>>>(n, c, h, x_dev, q_dev);
			}
			cudaDeviceSynchronize();
		}
			// cleanup and timings
		end = std::chrono::system_clock::now();
		cudaFree(x_dev);
		cudaFree(q_dev);
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << "GPU_source time (in ms): " << elapsed.count() << '\n';

			// eval
		start = std::chrono::system_clock::now();

		int nblocks_e = ceil(mdata/(float)nthreads); 
		dim3 threads_e(nthreads);
		dim3 blocks_e(nblocks_e);

		CUDA_CALL(cudaMalloc((void**) &y_dev, mdata*DIMENSIONS*sizeof(double)));
		CUDA_CALL(cudaMemcpy(y_dev, y, mdata*DIMENSIONS*sizeof(double), cudaMemcpyHostToDevice));
		cudaDeviceSynchronize();

		eval_ifgt_gpu<<<blocks_e,threads_e>>>(mdata, ncluster, h, y_dev, f_dev);
		cudaDeviceSynchronize();

		end = std::chrono::system_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << "GPU_target time (in ms): " << elapsed.count() << '\n';

		CUDA_CALL(cudaMemcpy(f, f_dev, mdata*sizeof(double), cudaMemcpyDeviceToHost));
		cudaDeviceSynchronize();

		cudaFree(y_dev);
		cudaFree(f_dev);

		hline();
		delete[] cidx;
		delete[] centers;
		delete[] np;
		delete[] clusters;
		
	}

	#endif	