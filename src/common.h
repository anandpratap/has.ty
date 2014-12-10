	#ifndef _COMMON_H
	#define _COMMON_H

	#define DEBUG 0
	#define INPUT_DATA_FILE "../matlab_test_codes/data.txt"
	#define OUTPUT_DATA_FILE "output.txt"

	
	#define PMAX 15
	#define NALPHA 120
	#define DIMENSIONS 2
	#define NMAX 131072
	#define KMAX 30

	struct cluster
	{
		double center[DIMENSIONS];
		double *x;
		double *q;
		unsigned int n;
		/* data */
	};

	
	#endif