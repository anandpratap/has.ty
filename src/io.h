	#ifndef _IO_H
	#define _IO_H
	// Contains input and output related functions

	/**
	* @brief Read in the data from the text file
	* Input file should have columns containing coordinates of data points followed by the exponential weights
	*
	* @param[in] ndata number of data points
	* @param[out] *x coordinates of the data points
	* @param[out] *q weight of the data points
	**/
	void read_data(int ndata, double *x, double *q);
	

	/**
	* @brief Write the final results into a file
	* output file contain columns containing exact value, approx value and absolute value of their difference
	*
	* @param[in] ndata number of data points
	* @param[in] *f_base benchmark gauss trasform 
	* @param[in] *f_approx approximate gauss transform
	**/
	
	void write_data(int ndata, double *f_base, double *f_approx);
	
	/**
	* @brief Parse the commandline inputs
	**/
	void parse_inputs(int argc, char **argv, int *ndata, int *mdata, int *nthreads, int *ncluster, bool *if_exact_gt, bool *if_host_ifgt, bool *if_device_ifgt);

	#endif