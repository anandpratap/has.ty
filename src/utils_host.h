	#ifndef _UTILS_H
	#define _UTILS_H
	
	#include "string.h"
	#include "common.h"
	
	
	/**
	* @brief Set first n elements of array x to zero
	* @param[in] n
	* @param[inout] *x
	**/
	void set_zeros(int n, double *x);

	/** 
	* @brief Return relative error as defined in the reference sum(abs(f_exact - f_approx))/Q
	* 
	* @param[in] ndata number of data points
	* @param[in] *f_base benchmark value
	* @param[in] *f_approx approximated value
	* @param[in] *q the weights of data points
	**/
	double calc_error(int ndata, double *f_base, double *f_approx, double *q);

	/**
	* @brief Return square of the L2 norm
	*
	* @param[in] *x 
	**/
	double l2normsq(double *x);

	unsigned int get_max_num_points(int ncluster, struct cluster *clusters);

	/**
	* @brief Draw a red color line
	**/
	void hline();

	/**
	* @brief printf in green color
	**/
	void gprintf(const char *x);
	
	#endif