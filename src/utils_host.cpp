	#ifndef _UTILS_H
	#define _UTILS_H
	
	#include <cmath> // exp
	#include <cassert>
	#include <fstream>
	#include "string.h"
	// local includes
	#include "common.h"
	#include "utils_host.h"
	#include <stdlib.h>
	
	void set_zeros(int n, double *x){
		memset(x, 0, n*sizeof(*x));
	}

	double calc_error(int ndata, double *f_base, double *f_approx, double *q){
		double error = 0.0;
		double Q = 0.0;
		for(int i=0; i < ndata; i++){
			error += fabs(f_base[i] - f_approx[i]);
			Q += fabs(f_base[i]);
		}
		error /= Q;
		return error;
	}

	

	double l2normsq(double *x){
		double norm = 0.0;
		for(int d=0; d<DIMENSIONS; d++){
			norm += x[d]*x[d];
		}
		return norm;
	}


	unsigned int get_max_num_points(int ncluster, struct cluster *clusters){
		unsigned int maxnp = 0;
		for(int c=0; c<ncluster; c++){
			if(clusters[c].n > maxnp){
				maxnp = clusters[c].n;
			}
		}
		return maxnp;
	}

	
	void hline(){
		printf("\x1b[31m---------------------------------\x1b[0m\n");
	}

	void gprintf(const char *x){
		printf("\x1b[32m%s\x1b[0m", x);
	}


	#endif