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
	__device__ void calc_monomial_gpu(int d, double *dx, int pmax, double *source_monomial){
		int *heads = new int[d]();
		int order = pmax - 1;

		source_monomial[0] = 1.0;
		
		int m, tail, head;
		m = 1;
		tail = 0;
		
		for(int ord=0; ord < order; ord++){
			for(int i=0; i < d; i++){

				head = heads[i];
				heads[i] = m;
				for(int j=head; j<=tail; j++){
					source_monomial[m] = dx[i]*source_monomial[j];
					m++;
				}

			}
			tail = m - 1;
		}

		delete[] heads;
	}


	// calculate coefficients
	__device__ void calc_constant_gpu(int d, int pmax, double *constant){
		int *heads = new int[d+1]();
		// THIS IS A TEMP FIX
		heads[d] = 1000000;

		int order = pmax - 1;
		constant[0] = 1.0;
		
		int *cinds = new int[nchoosek_gpu(pmax-1+d, d)]();
		cinds[0] = 0;

		int m, tail, head;
		m = 1;
		tail = 0;
		for(int ord=0; ord < order; ord++){
			for(int i=0; i < d; i++){
		
				head = heads[i];
				heads[i] = m;
				
				for(int j=head; j<=tail; j++){
					if(j < heads[i+1]){
						cinds[m] = cinds[j] + 1;
					}
					else{
						cinds[m] = 1;
					}
					constant[m] = 2*constant[j]/cinds[m];
					m++;
				}

			}
			tail = m - 1;
		}

		delete[] heads;
		delete[] cinds;
	}


	__device__ void calc_c_alpha_gpu(int n, int d, int pmax, double *x, double *c, double h, double *q, double *c_alpha){
		int nalpha = nchoosek_gpu(pmax-1+d, d);

		double *source_monomial = new double[nalpha]();
		double *constant = new double[nalpha]();
		double *dx = new double[d]();

		set_zeros_gpu(nalpha, c_alpha);
		
		for(int i = 0; i < n; i++){
			
			// calculate dx/h
			for(int j=0; j<d; j++){
				dx[j] = (x[i*d+j] - c[j])/h;
			}
			double dx2 = l2normsq_gpu(d, dx);

			// calc monomial contribution and add it to c_alpha
			calc_monomial_gpu(d, dx, pmax, source_monomial);
			for(int alpha=0; alpha < nalpha; alpha++){
				c_alpha[alpha] = c_alpha[alpha] + source_monomial[alpha]*q[i]*exp(-dx2);
			}

		}

		// calculate constant and multiple it with the monomials
		calc_constant_gpu(d, pmax, constant);
		for(int alpha=0; alpha < nalpha; alpha++){
			c_alpha[alpha] *= constant[alpha];
		}
		
		delete[] dx;
		delete[] constant;
		delete[] source_monomial;

	}


	__global__ void calc_ifgt_gauss_transform_gpu(int n, int m, int d, int pmax, double *x, double *y, double *c, double h, double *q, double *f){
		int nalpha = nchoosek_gpu(pmax-1+d, d);

		double *c_alpha = new double[nalpha]();
		double *target_monomial = new double[nalpha]();
		double *dy = new double[d]();
		
		// calcuate source coefficients
		calc_c_alpha_gpu(n, d, pmax, x, c, h, q, c_alpha);

		// for all test points
		for(int i=0; i< m; i++){
			
			set_zeros_gpu(nalpha, target_monomial);
			
			for(int j=0; j<d; j++){
				dy[j] = (y[i*d+j] - c[j])/h;
			}
			double dy2 = l2normsq_gpu(d, dy);
			
			// calculate target monomials
			calc_monomial_gpu(d, dy, pmax, target_monomial);

			// calculate the ifgt
			f[i] = 0.0;
			for(int alpha=0; alpha<nalpha; alpha++){
				f[i] += target_monomial[alpha]*c_alpha[alpha]*exp(-dy2);
			}

		}

		// clean up
		delete[] c_alpha;
		delete[] dy;
	}

	#endif