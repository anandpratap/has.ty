	#include <stdio.h>
	#include <stdlib.h>
	#include <iostream>
	#include <chrono>
	#include <string>
	#include "cmdline.h" // input parse
	
	#include "io.h"
	#include "utils.h"
	#include "ifgt.h"
	#include "common.h"
	

	#include "cuda.h"
	#include "cuda_runtime_api.h"
	#include "ifgt_gpu.h"
	
	
	
	// device constant memory variable
	__device__ __constant__ double constant_dev[NALPHA];
	__device__ __constant__ double centers_dev[KMAX*DIMENSIONS];

	int main(int argc, char **argv){

		int ndata, mdata, nthreads, ncluster;
		bool exact_gt, cpu_ifgt, gpu_ifgt;
		parse_inputs(argc, argv, &ndata, &mdata, &nthreads, &ncluster, &exact_gt, &cpu_ifgt, &gpu_ifgt);


		double h = 1.0;
		
		
	    // array data points
		// ith data point coordinates will have indices: i*d, i*d + 1, ..... i*d + d - 1;
		double *x = new double[ndata*DIMENSIONS];
		
		// coefficient of exponent
		double *q = new double[ndata];
		
		// read data
		read_data(ndata, x, q);
		
		double *f_exact = new double[mdata];

		if(exact_gt){

			hline();
			gprintf("CPU: Exact Gauss Transform\n");
			hline();
			auto start = std::chrono::system_clock::now();
			calc_exact_gauss_transform(ndata, mdata, DIMENSIONS, x, x, h, q, f_exact);	
			auto end = std::chrono::system_clock::now();
			auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			std::cout << "Exact time (in ms): " << elapsed.count() << '\n';
			hline();
			
		}

		// clustering ----------------
		double *centers = new double[DIMENSIONS*ncluster]();
		unsigned int *cidx = new unsigned int [ndata]();
		unsigned int *np = new unsigned int [ncluster]();
		struct cluster *clusters = new struct cluster[ncluster]();
		
		clustering(ndata, ncluster, x, q, cidx, centers, np, clusters);
		
		for(int i=0; i<ncluster; i++){
			printf("Cluster no: %d Number of points = %d\n", i, np[i]);
		}
		// ---------------------------
		

		double *f_approx_cpu = new double[mdata]();
		if(cpu_ifgt){
			hline();
			gprintf("CPU: Improved Fast Gauss Transform\n");
			hline();
			
			auto start = std::chrono::system_clock::now();
			calc_ifgt_gauss_transform(ndata, mdata, DIMENSIONS, ncluster, PMAX, x, x, centers, cidx, h, q, f_approx_cpu);
			auto end = std::chrono::system_clock::now();
			auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			std::cout << "CPU IFGT time (in ms): " << elapsed.count() << '\n';
		

			if(exact_gt){
		        // calculate and print error
				double error = calc_error(mdata, f_exact, f_approx_cpu, q);
				printf("CPU IFGT Error: %1.5E\n", error);
			}
		
			hline();
			
		}

		
		// write data to file f_exact f_approx abs(f_exact - fapprox)
		//write_data(n, f_exact, f_approx_cpu);


		
		if(gpu_ifgt){
		// cuda stuff
			hline();
			gprintf("GPU: Improved Fast Gauss Transform\n");
			hline();
		


			unsigned int max_num_points = get_max_num_points(ncluster, clusters);

			double *constant_cpu = new double[NALPHA]();
			calc_constant(constant_cpu);

			double *f_approx_cuda = new double[mdata]();
			double *x_dev;
			double *y_dev;
			double *q_dev;
			double *f_dev;
			
			
			printf("Max N point %d\n", max_num_points);
			GPU_ERROR_CHECK(cudaMalloc((void**) &x_dev, max_num_points*DIMENSIONS*sizeof(double)));
			GPU_ERROR_CHECK(cudaMalloc((void**) &centers_dev, ncluster*DIMENSIONS*sizeof(double)));
			GPU_ERROR_CHECK(cudaMalloc((void**) &q_dev, max_num_points*sizeof(double)));
			GPU_ERROR_CHECK(cudaMalloc((void**) &f_dev, mdata*sizeof(double)));
			GPU_ERROR_CHECK(cudaMemcpyToSymbol(constant_dev, constant_cpu, NALPHA*sizeof(double)));
			GPU_ERROR_CHECK(cudaMemcpyToSymbol(centers_dev, centers, ncluster*DIMENSIONS*sizeof(double)));

			
			auto start = std::chrono::system_clock::now();
				
			for(int c=0; c<ncluster; c++){
				struct cluster current_cluster = clusters[c];
				unsigned int n = current_cluster.n;

				GPU_ERROR_CHECK(cudaMemcpy(x_dev, current_cluster.x, n*DIMENSIONS*sizeof(double), cudaMemcpyHostToDevice));
				GPU_ERROR_CHECK(cudaMemcpy(q_dev, current_cluster.q, n*sizeof(double), cudaMemcpyHostToDevice));
				cudaDeviceSynchronize();

				int nblocks = ceil(current_cluster.n/(float)nthreads); 
				dim3 threads(nthreads);
				dim3 blocks(nblocks);

				calc_source_monomials_gpu<<<blocks,threads, nthreads*sizeof(double)>>>(n, c, h, x_dev, q_dev);
				cudaDeviceSynchronize();
			}
			// cleanup and timings
			auto end = std::chrono::system_clock::now();
			cudaFree(x_dev);
			cudaFree(q_dev);
			auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			std::cout << "GPU_source time (in ms): " << elapsed.count() << '\n';

			// eval
			start = std::chrono::system_clock::now();
			
			int nblocks = ceil(mdata/(float)nthreads); 
			dim3 threads(nthreads);
			dim3 blocks(nblocks);
			
			GPU_ERROR_CHECK(cudaMalloc((void**) &y_dev, mdata*DIMENSIONS*sizeof(double)));
			GPU_ERROR_CHECK(cudaMemcpy(y_dev, x, mdata*DIMENSIONS*sizeof(double), cudaMemcpyHostToDevice));
			cudaDeviceSynchronize();

			eval_ifgt_gpu<<<blocks,threads>>>(mdata, ncluster, h, y_dev, f_dev);
			cudaDeviceSynchronize();

			end = std::chrono::system_clock::now();
			elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			std::cout << "GPU_target time (in ms): " << elapsed.count() << '\n';

			GPU_ERROR_CHECK( cudaMemcpy(f_approx_cuda, f_dev, mdata*sizeof(double), cudaMemcpyDeviceToHost));
			cudaDeviceSynchronize();

			cudaFree(y_dev);
			cudaFree(centers_dev);
			cudaFree(f_dev);
	
		// calculate and print error
			if(exact_gt){
				double error = calc_error(mdata, f_exact, f_approx_cuda, q);
				printf("Error with respect to exact gauss transform %1.5E\n", error);
				write_data(mdata, f_exact, f_approx_cuda);
			}
			else if(cpu_ifgt){
				double error = calc_error(mdata, f_approx_cpu, f_approx_cuda, q);
				printf("Error with respect to cpu ifgt %1.5E\n", error);
				write_data(mdata, f_approx_cpu, f_approx_cuda);
			}

			else{
				gprintf("Nothing to compare against!\n");
			}
			
			delete[] f_approx_cuda;

			hline();
			
		}
		// clean up
		delete[] x;
		delete[] q;
		delete[] f_approx_cpu;
		delete[] f_exact;
		delete[] centers;
		delete[] clusters;
		delete[] np;
		delete[] cidx;

		return 0;
	}