	#ifndef _IFGT_HOST_H
	#define _IFGT_HOST_H
	#include <iostream>
	#include <chrono>
	#include <cmath>

	#include "common.h"
	#include "utils_host.h"
	#include "clustering_host.h"

	#include "ifgt_host.h"

		// calculate monomianls
	// used both for source and target
	void calc_monomial(int d, double *dx, int pmax, double *source_monomial){
		int heads[DIMENSIONS] = {0};
		
		source_monomial[0] = 1.0;
		
		int m, tail, head;
		m = 1;
		tail = 0;
		
		for(int ord=0; ord < PMAX - 1; ord++){
			for(int i=0; i < DIMENSIONS; i++){

				head = heads[i];
				heads[i] = m;
				for(int j=head; j<=tail; j++){
					source_monomial[m] = dx[i]*source_monomial[j];
					m++;
				}

			}
			tail = m - 1;
		}

	}


	// calculate coefficients
	void calc_constant(double *constant){
		int heads[DIMENSIONS+1] = {0};
		heads[DIMENSIONS] = std::numeric_limits<int>::max();

		constant[0] = 1.0;
		
		int cinds[NALPHA] = {0};
		cinds[0] = 0;

		int m, tail, head;
		m = 1;
		tail = 0;
		for(int ord=0; ord < PMAX-1; ord++){
			for(int i=0; i < DIMENSIONS; i++){
		
				head = heads[i];
				heads[i] = m;
				
				for(int j=head; j<=tail; j++){
					if(j < heads[i+1]){
						cinds[m] = cinds[j] + 1;
					}
					else{
						cinds[m] = 1;
					}
					constant[m] = 2*constant[j]/cinds[m];
					m++;
				}

			}
			tail = m - 1;
		}

		
	}


	void calc_c_alpha(int n, int d, int pmax, int nk, double *x, double *c, unsigned int *cidx, double h, double *q, double *c_alpha){


		double source_monomial[NALPHA] = {0};
		double constant[NALPHA] = {0};
		double dx[DIMENSIONS] = {0};

		set_zeros(NALPHA*nk, c_alpha);
		
		for(int i = 0; i < n; i++){
			
			// calculate dx/h
			for(int j=0; j<d; j++){
				dx[j] = (x[i*d+j] - c[cidx[i]*DIMENSIONS + j])/h;
			}
			double dx2 = l2normsq(dx);

			// calc monomial contribution and add it to c_alpha
			calc_monomial(d, dx, pmax, source_monomial);
			for(int alpha=0; alpha < NALPHA; alpha++){
				c_alpha[cidx[i]*NALPHA + alpha] += source_monomial[alpha]*q[i]*exp(-dx2);
			}

		}

		// calculate constant and multiple it with the monomials
		calc_constant(constant);
		for(int alpha=0; alpha < NALPHA; alpha++){
			for(int k=0; k<nk; k++){
				c_alpha[k*NALPHA + alpha] *= constant[alpha];
			}
		}
		
		
	}


	void calc_ifgt_gauss_transform(int n, int m, int d, int nk, int pmax, double *x, double *y, double *centers, unsigned int *cidx, double h, double *q, double *f){
		//double c_alpha[NALPHA] = {0};
		double *c_alpha = new double[NALPHA*nk]();
		double target_monomial[NALPHA] = {0};
		double dy[DIMENSIONS] = {0};
		
		// calcuate source coefficients
		calc_c_alpha(n, DIMENSIONS, pmax, nk, x, centers, cidx, h, q, c_alpha);
		printf("Source calculated\n");
		// for all test points
		for(int i=0; i< m; i++){
			// calculate the ifgt
			f[i] = 0.0;

			for(int k=0; k<nk; k++){
				set_zeros(NALPHA, target_monomial);
				for(int d=0; d<DIMENSIONS; d++){
					dy[d] = (y[i*DIMENSIONS+d] - centers[k*DIMENSIONS+d])/h;
				}

				double dy2 = l2normsq(dy);

			// calculate target monomials
				calc_monomial(DIMENSIONS, dy, pmax, target_monomial);
				for(int alpha=0; alpha<NALPHA; alpha++){
					f[i] += target_monomial[alpha]*c_alpha[k*NALPHA + alpha]*exp(-dy2);
				}
			}
		}
		delete[] c_alpha;
	}



	void ifgt_host(int ndata, int mdata, int ncluster, double h, double *x, double *q, double *y, double *f){
		hline();
		gprintf("CPU: Improved Fast Gauss Transform\n");
		hline();
		auto start = std::chrono::system_clock::now();

		// make clusters
		double *centers = new double[DIMENSIONS*ncluster]();
		unsigned int *cidx = new unsigned int [ndata]();
		unsigned int *np = new unsigned int [ncluster]();
		struct cluster *clusters = new struct cluster[ncluster]();
		
		clustering(ndata, ncluster, x, cidx, centers, np);
		rearrange_data(ndata, ncluster, x, q, cidx, centers, np, clusters);
		
		auto end = std::chrono::system_clock::now();
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << "CPU clustering time (in ms): " << elapsed.count() << '\n';



		start = std::chrono::system_clock::now();
		
		// calculate ifgt
		calc_ifgt_gauss_transform(ndata, mdata, DIMENSIONS, ncluster, PMAX, x, y, centers, cidx, h, q, f);
		
		end = std::chrono::system_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << "CPU IFGT time (in ms): " << elapsed.count() << '\n';
		hline();

		delete[] centers;
		delete[] cidx;
		delete[] np;
		delete[] clusters;
	}


	#endif