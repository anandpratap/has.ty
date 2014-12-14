	#ifndef _CLUSTERING_HOST_H
	#define _CLUSTERING_HOST_H

	#include <iostream>
	#include <fstream>
	#include <stdlib.h>
	#include <stdio.h>
	#include <cassert>
	
	#include "common.h"
	#include "utils_host.h"

	int diff(int n, unsigned int *cidx, unsigned int *cidx_old){
		int sum = 0;
		for(int i=0; i<n; i++){
			sum += abs(cidx[i] - cidx_old[i]);
		}
		return sum;
	}

	unsigned int get_minimum_idx(int nk, double *distance){
		unsigned int idx;
		double min_distance = 10000000000;
		for(int i=0; i<nk; i++){
			if(distance[i] < min_distance){
				min_distance = distance[i];
				idx = i;
			}
		}
		return idx;
	}
	
	void calc_centers(int n, int nk, double *x, unsigned int *cidx, double *centers, unsigned int *np){
		set_zeros(nk*DIMENSIONS, centers);
		for(int k=0; k < nk; k++){
			unsigned int npoints = 0;
			for(int i=0; i<n; i++){
				if(cidx[i] == k){
					for(int d=0; d<DIMENSIONS; d++){
						centers[k*DIMENSIONS + d] += x[DIMENSIONS*i + d];
					}	
					npoints += 1;
				}
			}
			np[k] = npoints;
			if(npoints != 0){
				for(int d=0; d<DIMENSIONS; d++){
					centers[k*DIMENSIONS + d] /= npoints;
				}
			}
		}
	}


	void clustering(int n, int nk, double *x, unsigned int *cidx, double *centers, unsigned int *np){
		// assign random cluster
		for(int i=0; i<n; i++){
			cidx[i] = rand() % nk;
		}

		double *distance = new double[nk]();
		double dx[DIMENSIONS];
		
		// calculate initial cluster center

		calc_centers(n, nk, x, cidx, centers, np);

		int iter = 0;
		unsigned int tmp_idx;
		while (1){
			// calculate distance and assign new clusters
			unsigned int delta = 0;
		
			for(int i=0; i< n; i++){
				for(int k=0; k<nk; k++){
					for(int d=0; d<DIMENSIONS; d++){
						dx[d] = x[i*DIMENSIONS+d] - centers[k*DIMENSIONS + d];
					}	
					distance[k] = l2normsq(dx);
				}
				tmp_idx = get_minimum_idx(nk, distance);
				if(tmp_idx != cidx[i]){
					delta += 1;
					cidx[i] = tmp_idx;
				}
				//printf("%d\n", cidx[i]);
			}

			// new cluster center
			calc_centers(n, nk, x, cidx, centers, np);

			iter += 1;
			if(delta == 0){
				printf("clustering iter: %d\n", iter);
				break;
			}
		}
		
		delete[] distance;
		
	}

	void rearrange_data(int ndata, int ncluster, double *x, double *q, unsigned int *cidx, double *centers, unsigned int *np, struct cluster *clusters){
		// assign clusters to struct
		for(int k=0; k<ncluster; k++){
			clusters[k].n = np[k];
			clusters[k].x = new double[clusters[k].n*DIMENSIONS]();
			clusters[k].q = new double[clusters[k].n]();
			for(int d=0; d<DIMENSIONS; d++){
				clusters[k].center[d] = centers[k*DIMENSIONS + d];
			}
			for(int i=0, j=0; i<ndata; i++){
				if(cidx[i] == k){
					for(int d=0; d<DIMENSIONS;d++){
						clusters[k].x[j*DIMENSIONS + d] = x[i*DIMENSIONS + d];
					}
					clusters[k].q[j] = q[i];
					j++;
				}
			}
		}
		
		std::fstream outfile("clustering.txt", std::ios_base::out);
		// this should be removed 
		for(int k=0; k<ncluster; k++){
			for(int i=0; i<clusters[k].n; i++){
				outfile << clusters[k].x[i*DIMENSIONS] << " " << clusters[k].x[i*DIMENSIONS+1] << " " <<clusters[k].center[0]<<" "<<clusters[k].center[1]<<" "<< k << "\n";
			}
		}

	}


	#endif