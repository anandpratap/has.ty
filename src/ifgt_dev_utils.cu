	#ifndef _IFGT_GPU_H
	#define _IFGT_GPU_H

	#include <stdlib.h>
	#include <stdio.h>
	#include <limits>  // max int
	#include <cmath> // exp

	// local includes	
	#include "utils_dev.h"
	#include "common.h" // DEBUG 
	
	// global device variables
	extern __device__ __constant__ double constant_dev[NALPHA];
	extern __device__ __constant__ double centers_dev[KMAX*DIMENSIONS];

	__device__ double coeff_global[KMAX*NALPHA] = {0};
	
	__device__ void AtomicAddDouble (volatile double *address, double value){
		
		unsigned long long oldval, newval, readback; 
		oldval = __double_as_longlong(*address);
		newval = __double_as_longlong(__longlong_as_double(oldval) + value);
		
		while ((readback=atomicCAS((unsigned long long *)address, oldval, newval)) != oldval)
		{
			oldval = readback;
			newval = __double_as_longlong(__longlong_as_double(oldval) + value);
		}
	}

	
	// calculate monomianls
	// used both for source and target
	__device__ inline void calc_monomial_gpu(int idx, int d, int pmax, int nalpha, double *dx, double *source_monomial){
		unsigned short int heads[DIMENSIONS] = {0};

		source_monomial[0] = 1.0;
		unsigned short int m, tail, head;
		m = 1;
		tail = 0;
		for(unsigned char ord=0; ord < PMAX - 1; ord++){
			for(unsigned char i=0; i < DIMENSIONS; i++){
				head = heads[i];
				heads[i] = m;
				for(unsigned short int j=head; j<=tail; j++){
					source_monomial[m] = dx[i]*source_monomial[j];
					m++;
				}

			}
			tail = m - 1;
		}

	}

	


	__global__ void calc_source_monomials_gpu(int n, int c, double h, double *x, double *q){
		// get ids and indices
		unsigned short int tid = threadIdx.x;
		unsigned short int blockdim = blockDim.x;
		unsigned int idx = blockIdx.x*blockdim + tid;
		extern __shared__ double coeff_alpha[];
		
		__syncthreads();
		double dx[DIMENSIONS] = {0};
		double ex;
		double source_monomial_local[NALPHA] = {0};

		// if index in bounds
		if(idx < n){
			
			// calculate dx/h
			for(unsigned short int d=0; d<DIMENSIONS; d++){
				dx[d] = (x[idx*DIMENSIONS+d] - centers_dev[c*DIMENSIONS + d])/h;
			}
		    // calculate source monomial
			calc_monomial_gpu(idx, DIMENSIONS, PMAX, NALPHA, dx, source_monomial_local);
			ex = exp(-l2normsq_dev(dx))*q[idx];

			// source*exp(-dx^2/h^2),  2^{alpha}/alpha! is multiplied in the reduction step to save the computations
			for(unsigned short int alpha=0; alpha<NALPHA; alpha++){
				source_monomial_local[alpha] = source_monomial_local[alpha]*ex*constant_dev[alpha];
			}
		}
		__syncthreads();
		// reduce here 
		for(unsigned short int alpha=0; alpha<NALPHA; alpha++){
			coeff_alpha[tid] = source_monomial_local[alpha];
			__syncthreads();
			

			for(unsigned short int i=blockdim/2; i>=1; i >>= 1){
				if(tid < i){
					coeff_alpha[tid] += coeff_alpha[tid + i];
				}
				__syncthreads();
			}

			// does not work without debug flag!!
			/*if(tid < 32){
				coeff_alpha[tid] += coeff_alpha[tid + 32];
				coeff_alpha[tid] += coeff_alpha[tid + 16];
				coeff_alpha[tid] += coeff_alpha[tid + 8];
				coeff_alpha[tid] += coeff_alpha[tid + 4];
				coeff_alpha[tid] += coeff_alpha[tid + 2];
				coeff_alpha[tid] += coeff_alpha[tid + 1];
			}*/

			__syncthreads();

			if(tid == 0 ){
				AtomicAddDouble(&coeff_global[c*NALPHA + alpha], coeff_alpha[0]);
			}
			__syncthreads();
		}
		__syncthreads();
	}

	__global__ void eval_ifgt_gpu(int mdata, int ncluster, double h, double *y, double *f){
		// get indices
		unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x; 

		// local variables
		double dy[DIMENSIONS] = {0};
		double target_monomial_local[NALPHA] = {0};
		
		// using tmp so that we dont write to global memory everytime
		double tmp = 0.0;
		double ex;

		if(idx < mdata){
			
			for(unsigned short int c=0;c<ncluster;c++){
				// calculate dy/h
				for(int d=0; d<DIMENSIONS; d++){
					dy[d] = (y[idx*DIMENSIONS + d] - centers_dev[c*DIMENSIONS + d])/h;
				}
				ex = exp(-l2normsq_dev(dy));
		    // calc target monomial
				calc_monomial_gpu(idx, DIMENSIONS, PMAX, NALPHA, dy, target_monomial_local);
				for(unsigned short int alpha=0; alpha<NALPHA; alpha++){
					tmp += target_monomial_local[alpha]*coeff_global[c*NALPHA+alpha]*ex;
				}
			}
			f[idx] = tmp;
		}
	}


	#endif