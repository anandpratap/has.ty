	#ifndef _IFGT_H
	#define _IFGT_H

	// Calculate the monomial given dx, note that the input dx include division by h i.e. actual dx/h
	// It can be used to calculate both the target and source monomial
	// d -> number of dimensions
	// dx -> x - c where x is a data point and c is cluster center
	//
	// source_mononial -> calculate monomials
	
	void calc_monomial(int d, double *dx, int pmax, double *source_monomial);

	// Calculate 2^alpha/(alpha!)
	// d -> number of dimensions
	// pmax -> maximum number of terms
	// 
	// constant -> output coefficients of monomials
	void calc_constant(double *constant);

	// Calculate C_alpha
	// x -> number of data points
	// d -> number of dimensions
	// pmax -> maximum number of terms
	// x -> data points
	// c -> cluster centers
	// h -> bandwidth
	// q -> weights
	//
	// c_alpha -> calculated coefficients
	void calc_c_alpha(int n, int d, int pmax, double *x, double *c, double h, double *q, double *c_alpha);
	
	// Calculate ifgt gauss tranform form m test points given by y
	void calc_ifgt_gauss_transform(int n, int m, int d, int nk, int pmax, double *x, double *y, double *centers, unsigned int *cidx, double h, double *q, double *f);

	#endif