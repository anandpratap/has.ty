	#ifndef _IO_H
	#define _IO_H
	// Contains input and output related functions

	
	// read the data from the text file
	// input file should have columns containing coordinates of data points followed by the exponential weights
	// n -> number of data points
	// d -> number of dimensions
	// x -> coodinates
	// q -> weights
	void read_data(int n, int d, double *x, double *q);
	

	// write the final results into a file
	// output file contain columns containing exact value, approx value and absolute value of their difference
	void write_data(int n, double *f_exact, double *f_approx);
	
	#endif