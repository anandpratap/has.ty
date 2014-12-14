	#ifndef _CLUSTERING_HOST_H
	#define _CLUSTERING_HOST_H

	#include <iostream>
	#include <stdlib.h>
	#include <stdio.h>
	
	/**
	* @brief Calculate cluster centers
	* 
	* @param[in] ndata number of data points
	* @param[in] ncluster number of clusters
	* @param[in] *x contains the coordinates of data points
	* @param[in] *cidx array containing cluster number for every data points
	* @param[out] *centers pointer to array containing centers
	* @param[out] *np array containing number of data points for all the clusters 
	**/
	void calc_centers(int ndata, int ncluster, double *x, unsigned int *cidx, double *centers, unsigned int *np);
	

	/**
	* @brief Perform k-mean clustering on host
	*
	* @param[in] ndata number of data points
	* @param[in] ncluster number of clusters
	* @param[in] *x contains the coordinate of data points
	* @param[out] *cidx array containing the cluster id
	* @param[out] *centers contains cluster centers
	* @param[out] *np containing number of points in a cluster
	**/
	void clustering(int ndata, int ncluster, double *x, unsigned int *cidx, double *centers, unsigned int *np);

	/**
	* @brief: Rearrange data in the form of struct so that it is fed to the device in memory contiguous form
	* 
	* @param[in] ndata number of data points
	* @param[in] ncluster number of clusters
	* @param[in] *x contains the coordinate of data points
	* @param[in] *q contains the weight of data points
	* @param[in] *cidx array containing the cluster id
	* @param[in] *centers contains cluster centers
	* @param[in] *np containing number of points in a cluster
	* @param[out] *clusters array of struct of clusters, each cluster contains its center, pointer to data points and number of points
	**/
	void rearrange_data(int ndata, int ncluster, double *x, double *q, unsigned int *cidx, double *centers, unsigned int *np, struct cluster *clusters);
	
	#endif