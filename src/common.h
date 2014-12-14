	#ifndef _COMMON_H
	#define _COMMON_H
	#include <iostream>
	#define DEBUG 1
	#define INPUT_DATA_FILE "../matlab_test_codes/data.txt"
	#define OUTPUT_DATA_FILE "output.txt"

	#define MAX_DOUBLE 1e100
	#define PMAX 15
	#define NALPHA 120
	#define DIMENSIONS 2
	#define NMAX 131072
	#define KMAX 20

	struct cluster{
		double center[DIMENSIONS];
		double *x;
		double *q;
		unsigned int n;
	};

	#endif