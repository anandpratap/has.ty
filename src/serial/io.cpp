	#ifndef _IO_H
	#define _IO_H
	
	#include <stdio.h>
	#include <fstream>
	#include <cmath>
	// local includes
	#include "common.h" // DEBUG
	
	
	void read_data(int n, int d, double *x, double *q){
		std::fstream inpfile(INPUT_DATA_FILE, std::ios_base::in);
		inpfile.precision(16);
		
		for(int i=0; i<n; i++){
			for(int j=0; j< d; j++){
				inpfile >> x[d*i + j];
			}
			inpfile >> q[i];
		}
	
	}

	void write_data(int n, double *f_exact, double *f_approx){
		std::fstream outfile(OUTPUT_DATA_FILE, std::ios_base::out);
		// set precision close to double's
		outfile.precision(16);
		
		for(int i=0; i < n; i++){
			outfile << f_exact[i] << " " << f_approx[i] << " " << fabs(f_exact[i]-f_approx[i]) << "\n";
		}

	}

	#endif