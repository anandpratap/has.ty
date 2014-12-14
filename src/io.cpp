	#ifndef _IO_H
	#define _IO_H
	
	#include <stdio.h>
	#include <fstream>
	#include <cmath>
	// local includes
	#include "common.h" // DEBUG
	#include "cmdline.h"
	
	void read_data(int ndata, double *x, double *q){
		std::fstream inpfile(INPUT_DATA_FILE, std::ios_base::in);
		inpfile.precision(16);
		double tmp;

		for(int i=0; i<ndata; i++){
			for(int d=0; d< DIMENSIONS; d++){
				inpfile >> x[DIMENSIONS*i + d];
			}
			inpfile >> q[i];
			inpfile >> tmp;
		}
	}

	void write_data(int ndata, double *f_base, double *f_approx){
		std::fstream outfile(OUTPUT_DATA_FILE, std::ios_base::out);
		// set precision close to double's
		outfile.precision(16);
		
		for(int i=0; i < ndata; i++){
			outfile << f_base[i] << " " << f_approx[i] << " " << fabs(f_base[i]-f_approx[i]) << "\n";
		}

	}


	void parse_inputs(int argc, char **argv, int *ndata, int *mdata, int *nthreads, int *ncluster, bool *if_exact_gt, bool *if_host_ifgt, bool *if_device_ifgt){
		
		cmdline::parser a;
		
		a.add<int>("ndata", 'd', "number of data points", true, 10000);
		a.add<int>("mdata", 'm', "number of test points", true, 10000);
		a.add<int>("nthreads", 't', "number of threads per block", false, 128);
		a.add<int>("ncluster", 'c', "number of clusters", false, 1, cmdline::range(1, 100));
		a.add("exact_gt", '\0', "calculate exact gauss transform");
		a.add("host_ifgt", '\0', "calculate improved fast gauss transform on cpu");
		a.add("device_ifgt", '\0', "calculate improved fast gauss transform on gpu");
		a.parse_check(argc, argv);
		
		
		*ndata = a.get<int>("ndata");
		*mdata = a.get<int>("mdata");
		*nthreads = a.get<int>("nthreads");
		*ncluster = a.get<int>("ncluster");
		
		*if_exact_gt = (a.exist("exact_gt")) ? true : false;
		*if_host_ifgt = (a.exist("host_ifgt")) ? true : false;
		*if_device_ifgt = (a.exist("device_ifgt")) ? true : false;
		
	}

	#endif