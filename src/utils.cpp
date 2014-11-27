	#ifndef _UTILS_H
	#define _UTILS_H
	
	#include <cmath> // exp
	#include "string.h"
	// local includes
	#include "common.h"
	#include "utils.h"
	
	template <class dtype>
	void set_zeros(int n, dtype *x){
		memset(x, 0, n*sizeof(*x));
	}


	int nchoosek(int n, int k){
		int n_k = n - k;
		if (k < n_k){
			k = n_k;
			n_k = n - k;
		}
		int  nchsk = 1; 
		for ( int i = 1; i <= n_k; i++){
			nchsk *= (++k);
			nchsk /= i;
		}

		return nchsk;
	}

	void calc_centroid(int n, int d, double *x, double *c){
		set_zeros(d, c);

		for(int i=0; i<n; i++){
			for(int j=0; j<d; j++){
				c[j] += x[d*i +j];
			}
		}

		for(int i=0; i < d; i++){
			c[i] = c[i]/n;
		}
	}

	double calc_error(int n, double *f, double *f_approx, double *q){
		double error = 0.0;
		double Q = 0.0;
		
		for(int i=0; i < n; i++){
			error += fabs(f[i] - f_approx[i]);
			Q += fabs(q[i]);
		}

		error /= Q;

		return error;
	}

	void calc_exact_gauss_transform(int n, int m, int d, double *x, double *y, double h, double *q, double *f){
		set_zeros(m, f);
		for(int i=0; i < m; i++){
			for(int j=0; j < n; j++){
				double dx2 = pow(x[j*d] - y[i*d], 2) + pow(x[j*d+1] - y[i*d+1], 2);
				dx2 /= h*h;
				f[i] += q[j]*exp(-dx2);
			}
		}
	}


	double l2normsq(int n, double *x){
		double norm = 0.0;
		for(int i=0; i<n; i++){
			norm += x[i]*x[i];
		}
		return norm;
	}

	#endif