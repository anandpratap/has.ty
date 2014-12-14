	#ifndef _EXACT_GT_H
	#define _EXACT_GT_H
	#include <cmath>
	#include <stdlib.h>

	#include "exact_gt.h"
	#include "utils_host.h"
	#include "common.h"

	void exact_gt(int ndata, int mdata, double h, double *x, double *q, double *y, double *f){
		set_zeros(mdata, f);
		for(int i=0; i < mdata; i++){
			for(int j=0; j < ndata; j++){
				double dx2 = 0.0;
				for(int d=0; d<DIMENSIONS; d++){
					dx2 += pow(x[j*DIMENSIONS + d] - y[i*DIMENSIONS + d], 2);
				}
				dx2 /= h*h;
				f[i] += q[j]*exp(-dx2);
			}
		}
	}

	#endif