	#ifndef _UTILS_H
	#define _UTILS_H
	
	#include "string.h"

	// set first n elements of array x to zero
	template <class dtype>
	void set_zeros(int n, dtype *x);

	// calculate n!/k!(n-k)!
	int nchoosek(int n, int k);

	// calculate centroi of given points
	void calc_centroid(int n, int d, double *x, double *c);

	// calculate relative error as defined in the reference sum(abs(f_exact - f_approx))/Q
	double calc_error(int n, double *f, double *f_approx, double *q);

	// calculate exact gauss tranform
	void calc_exact_gauss_transform(int n, int m, int d, double *x, double *y, double h, double *q, double *f);

	// calculate square of the l2 norm
	double l2normsq(int n, double *x);

	#endif