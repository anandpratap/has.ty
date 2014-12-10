	#ifndef _IO_H
	#define _IO_H
	// Contains input and output related functions

	
	// read the data from the text file
	// input file should have columns containing coordinates of data points followed by the exponential weights
	// n -> number of data points
	// d -> number of dimensions
	// x -> coodinates
	// q -> weights
	void read_data(int n, double *x, double *q);
	

	// write the final results into a file
	// output file contain columns containing exact value, approx value and absolute value of their difference
	void write_data(int n, double *f_exact, double *f_approx);
	

	void parse_inputs(int argc, char **argv, int *ndata, int *mdata, int *nthreads, int *ncluster, bool *exact_gt, bool *cpu_ifgt, bool *gpu_ifgt);

	#endif