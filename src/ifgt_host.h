	#ifndef _IFGT_HOST_H
	#define _IFGT_HOST_H
	#include <iostream>
	#include "clustering_host.h"
	
	/**
	* @brief Calculate the monomial coefficients
	*
	* @param[out] *constant array containg the monomial coefficients
	**/
	void calc_constant(double *constant);

	/**
	* @brief Calculate the improved fast gauss transform on the host
	* 
	* @param[in] ndata number of data points
	* @param[in] mdata number of test points
	* @param[in] ncluster number of clusters
	* @param[in] h bandwidth
	* @param[in] *x contains the coordinate of data points
	* @param[in] *q contains the weight of data points
	* @param[in] *y contains the coordinate of test points
	* @param[out] *f contains the gauss tranform at test points
	**/
	void ifgt_host(int ndata, int mdata, int ncluster, double h, double *x, double *q, double *y, double *f);
	

	#endif