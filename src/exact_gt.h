	#ifndef _EXACT_GT_H
	#define _EXACT_GT_H
	
	#include "exact_gt.h"
	#include "utils_host.h"
	
	/**
	* @brief Calculate the exact gauss transform on the host
	* 
	* @param[in] ndata number of data points
	* @param[in] mdata number of test points
	* @param[in] h bandwidth
	* @param[in] *x contains the coordinate of data points
	* @param[in] *q contains the weight of data points
	* @param[in] *y contains the coordinate of test points
	* @param[out] *f contains the gauss tranform at test points
	**/
	
	void exact_gt(int ndata, int mdata, double h, double *x, double *q, double *y, double *f);

	#endif