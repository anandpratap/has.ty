	#include <stdio.h>
	#include <stdlib.h>
	#include <iostream>
	#include <chrono>
	#include "io.h"
	#include "utils.h"
	//#include "utils_gpu.cu"
	#include "ifgt.h"
	#include "common.h"
	#include "cuda.h"
	#include "cuda_runtime_api.h"
	#include "ifgt_gpu.h"
	int main(void){
		// ifgt parameters
		// n <int> -> number of data points
		// d <int> -> number of dimensions
		// pmax <int> -> max number of terms in the expansion
		// h <double> -> bandwidth


		int n = 10000;
		int d = 2;
		int pmax = 25;
		double h = 1.0;
		
		// number of test points, for this case this is same as data points
		int m = n;


		// array data points
		// ith data point coordinates will have indices: i*d, i*d + 1, ..... i*d + d - 1;
		double *x = new double[n*d];
		
		// coefficient of exponent
		double *q = new double[n];
		
		// read data
		read_data(n, d, x, q);

		
		// calculate exact gauss transform
		double *f_exact = new double[m];
		auto start = std::chrono::system_clock::now();
		calc_exact_gauss_transform(n, n, d, x, x, h, q, f_exact);
		auto end = std::chrono::system_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << "Exact time (in ms): " << elapsed.count() << '\n';
		

		// cluster center
		double *c = new double[d]();
		calc_centroid(n, d, x, c);
		
		// calculate ifgt
		double *f_approx = new double[m];
		start = std::chrono::system_clock::now();
		calc_ifgt_gauss_transform(n, n, d, pmax, x, x, c, h, q, f_approx);
		end = std::chrono::system_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << "IFGT time (in ms): " << elapsed.count() << '\n';
		
		// calculate and print error
		double error = calc_error(n, f_exact, f_approx, q);
		printf("%1.5E\n", error);
		// write data to file f_exact f_approx abs(f_exact - fapprox)
		write_data(n, f_exact, f_approx);


		// cuda stuff
		double *x_dev;
		double *c_dev;
		double *q_dev;
		double *f_approx_dev;
		cudaMalloc((void**) &x_dev, n*d*sizeof(double));
		cudaMalloc((void**) &c_dev, d*sizeof(double));
		cudaMalloc((void**) &q_dev, n*sizeof(double));
		cudaMalloc((void**) &f_approx_dev, m*sizeof(double));
		cudaMemcpy(x_dev, x, n*d*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(c_dev, c, d*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(q_dev, q, n*sizeof(double), cudaMemcpyHostToDevice);
		start = std::chrono::system_clock::now();
		calc_ifgt_gauss_transform_gpu<<<1,1>>>(n, n, d, pmax, x_dev, x_dev, c_dev, h, q_dev, f_approx_dev);
		cudaDeviceSynchronize();
		end = std::chrono::system_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << "GPU time (in ms): " << elapsed.count() << '\n';
		cudaMemcpy(f_approx, f_approx_dev, m*sizeof(double), cudaMemcpyDeviceToHost);
		cudaFree(x_dev);
		cudaFree(c_dev);
		cudaFree(q_dev);
		cudaFree(f_approx_dev);

		// calculate and print error
		error = calc_error(n, f_exact, f_approx, q);
		printf("%1.5E\n", error);
		write_data(n, f_exact, f_approx);
		// clean up
		delete[] x;
		delete[] q;
		delete[] f_exact;
		delete[] f_approx;
		return 0;

	}