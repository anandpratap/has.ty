	#ifndef _IFGT_DEV_H
	#define _IFGT_DEV_H
	
	/**
	* @brief Calculate the improved fast gauss transform on the device
	* 
	* @param[in] ndata number of data points
	* @param[in] mdata number of test points
	* @param[in] ncluster number of clusters
	* @param[in] h bandwidth
	* @param[in] *x contains the coordinate of data points
	* @param[in] *q contains the weight of data points
	* @param[in] *y contains the coordinate of test points
	* @param[out] *f contains the gauss tranform at test points
	* @param[in] nthreads number of threads
	**/
	void ifgt_dev(int ndata, int mdata, int ncluster, double h, double *x, double *q, double *y, double *f, int nthreads);

	#endif