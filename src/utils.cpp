	#ifndef _UTILS_H
	#define _UTILS_H
	
	#include <cmath> // exp
	#include <fstream>
	#include "string.h"
	// local includes
	#include "common.h"
	#include "utils.h"
	#include <stdlib.h>
	
	void set_zeros(int n, double *x){
		memset(x, 0, n*sizeof(*x));
	}


	int nchoosek(int n, int k){
		int n_k = n - k;
		if (k < n_k){
			k = n_k;
			n_k = n - k;
		}
		int  nchsk = 1; 
		for ( int i = 1; i <= n_k; i++){
			nchsk *= (++k);
			nchsk /= i;
		}

		return nchsk;
	}

	void calc_centroid(int n, int d, double *x, double *c){
		set_zeros(d, c);

		for(int i=0; i<n; i++){
			for(int j=0; j<d; j++){
				c[j] += x[d*i +j];
			}
		}

		for(int i=0; i < d; i++){
			c[i] = c[i]/n;
		}
	}

	double calc_error(int n, double *f, double *f_approx, double *q){
		double error = 0.0;
		double Q = 0.0;
		
		for(int i=0; i < n; i++){
			error += fabs(f[i] - f_approx[i]);
			Q += fabs(f[i]);
		}

		error /= Q;

		return error;
	}

	void calc_exact_gauss_transform(int n, int m, int d, double *x, double *y, double h, double *q, double *f){
		set_zeros(m, f);
		for(int i=0; i < m; i++){
			for(int j=0; j < n; j++){
				double dx2 = pow(x[j*d] - y[i*d], 2) + pow(x[j*d+1] - y[i*d+1], 2);
				dx2 /= h*h;
				f[i] += q[j]*exp(-dx2);
			}
		}
	}


	double l2normsq(int n, double *x){
		double norm = 0.0;
		for(int i=0; i<n; i++){
			norm += x[i]*x[i];
		}
		return norm;
	}

	int diff(int n, unsigned int *cidx, unsigned int *cidx_old){
		int sum = 0;
		for(int i=0; i<n; i++){
			sum += (cidx[i] - cidx_old[i]);
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


	void calc_centers(int n, int nk, double *centers, unsigned int *cidx, double *x, unsigned int *np){
		unsigned int npoints;
		for(int k=0; k < nk; k++){
			for(int d=0; d<DIMENSIONS; d++){
				centers[k*DIMENSIONS + d] = 0.0;
			}
			npoints = 0;
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


	void clustering(int n, int nk, double *x, double *q, unsigned int *cidx, double *centers, unsigned int *np, struct cluster *clusters){
		// assign random cluster
		for(int i=0; i<n; i++){
			cidx[i] = rand() % nk;
		}
		
		printf("%p\n", cidx);
		printf("%p\n", centers);

		unsigned int *cidx_old = new unsigned int[n]();
		double *distance = new double[nk]();
		double dx[DIMENSIONS];
		unsigned int npoints;
		// calculate initial cluster center

		calc_centers(n, nk, centers, cidx, x, np);

		

		int iter = 0;
		while (1){
			for(int i=0; i<n; i++){
				cidx_old[i] = cidx[i];
			}
			// calculate distance and assign new clusters
			for(int i=0; i< n; i++){
				for(int k=0; k<nk; k++){
					for(int d=0; d<DIMENSIONS; d++){
						dx[d] = x[i*DIMENSIONS+d] - centers[k*DIMENSIONS + d];
					}	
					distance[k] = l2normsq(DIMENSIONS, dx);
				}
				cidx[i] = get_minimum_idx(nk, distance);
				//printf("%d\n", cidx[i]);
			}

			// new cluster center
			calc_centers(n, nk, centers, cidx, x, np);
			iter += 1;
			if(diff(n, cidx, cidx_old) == 0){
				printf("clustering iter: %d\n", iter);
				break;
			}
		}

		delete[] distance;
		std::fstream outfile("clustering.txt", std::ios_base::out);
				
	/*
	*/	

		for(int k=0; k<nk; k++){
			clusters[k].n = np[k];
			clusters[k].x = new double[clusters[k].n*DIMENSIONS]();
			clusters[k].q = new double[clusters[k].n]();
			for(int d=0; d<DIMENSIONS; d++){
				clusters[k].center[d] = centers[k*DIMENSIONS + d];
			}
			for(int i=0, j=0; i<n; i++){
				if(cidx[i] == k){
					for(int d=0; d<DIMENSIONS;d++){
						clusters[k].x[j*DIMENSIONS + d] = x[i*DIMENSIONS + d];
					}
					clusters[k].q[j] = q[i];
					j++;
				}
			}
		}
		
		for(int k=0; k<nk; k++){
			for(int i=0; i<clusters[k].n; i++){
				outfile << clusters[k].x[i*DIMENSIONS] << " " << clusters[k].x[i*DIMENSIONS+1] << " " <<clusters[k].center[0]<<" "<<clusters[k].center[1]<<" "<< k << "\n";

			}
		}
		

	}


	unsigned int get_max_num_points(int ncluster, struct cluster *clusters){
		unsigned int maxnp = 0;
		for(int c=0; c<ncluster; c++){
			if(clusters[c].n > maxnp){
				maxnp = clusters[c].n;
			}
		}
		return maxnp;
	}

	
void hline(){
		printf("\x1b[31m---------------------------------\x1b[0m\n");
	}

	void gprintf(const char *x){
		printf("\x1b[32m%s\x1b[0m", x);
	}


	#endif