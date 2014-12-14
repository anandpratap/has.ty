	#include <stdio.h>
	#include <stdlib.h>
	#include <iostream>
	#include <fstream>
	#include <chrono>
	#include <string>
	#include "cmdline.h" // input parse
	#include <cassert>

	#include "io.h"
	#include "utils_host.h"
	#include "common.h"
	

	
	#include "ifgt_host.h"
	#include "ifgt_dev.h"
	#include "exact_gt.h"

	int main(int argc, char **argv){

		int ndata, mdata, nthreads, ncluster;
		bool if_exact_gt, if_host_ifgt, if_dev_ifgt;
		parse_inputs(argc, argv, &ndata, &mdata, &nthreads, &ncluster, &if_exact_gt, &if_host_ifgt, &if_dev_ifgt);


		// read in inputs and stuff
		double h = 1.0;
		
		// allocate and read in data points & coeffs
		double *x = new double[ndata*DIMENSIONS];
		double *q = new double[ndata];
		read_data(ndata, x, q);
		
		// read in the test points, for now then are same as data points
		double *y = new double[ndata*DIMENSIONS];
		for(int i=0; i<ndata*DIMENSIONS; i++)
			y[i] = x[i];


		//
		// EXACT GAUSS TRANSFORM ~~~~~
		//
		double *f_exact = new double[mdata];
		
		if(if_exact_gt){
			hline();
			gprintf("CPU: Exact Gauss Transform\n");
			hline();
			auto start = std::chrono::system_clock::now();
			exact_gt(ndata, mdata, h, x, q, y, f_exact);	
			auto end = std::chrono::system_clock::now();
			auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			std::cout << "Exact time (in ms): " << elapsed.count() << '\n';
			hline();
			
		}
		//
		// EXACT GAUSS TRANSFORM ENDS ~~~~~
		//

		
		//
		// HOST IMPROVED FAST GAUSS TRANSFORM ~~~~~
		//
		double *f_approx_host = new double[mdata]();

		if(if_host_ifgt){
			ifgt_host(ndata, mdata, ncluster, h, x, q, y, f_approx_host);
			if(if_exact_gt){
				double error = calc_error(mdata, f_exact, f_approx_host, q);
				printf("Host IFGT Error: %1.5E\n", error);
			}
		}
		//
		// HOST IMPROVED FAST GAUSS TRANSFORM ENDS ~~~~~
		//

		//
		// DEVICE IMPROVED FAST GAUSS TRANSFORM ~~~~~
		//
		double *f_approx_dev = new double[mdata]();
		
		if(if_dev_ifgt){
			ifgt_dev(ndata, mdata, ncluster, h, x, q, y, f_approx_dev, nthreads);
			
			if(if_exact_gt){
				double error = calc_error(mdata, f_exact, f_approx_dev, q);
				printf("Device IFGT Error wrt exact: %1.5E\n", error);
			}
			else if(if_host_ifgt){
				double error = calc_error(mdata, f_approx_host, f_approx_dev, q);
				printf("Device IFGT Error wrt host: %1.5E\n", error);

			}

		}
	    //
		// DEVICE IMPROVED FAST GAUSS TRANSFORM ENDS ~~~~~
		//
		
		// clean up
		delete[] x;
		delete[] q;
		delete[] y;
		delete[] f_exact;
		delete[] f_approx_host;
		delete[] f_approx_dev;

		return 0;
	}