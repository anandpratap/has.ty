__global__ void calc_source_monomials_gpu(int n, int nk, int pmax, int nalpha, double h, double *x, double *centers, unsigned int *cidx, double *q, double *c_alpha_dev){
		// get ids and indices
		int tid = threadIdx.x;
		int blockdim = blockDim.x;
		int idx = blockIdx.x*blockdim + tid;
		__shared__ double c_ak[KMAX*NALPHA];
		double dx2;
		double dx[DIMENSIONS] = {0};
		double ex;
		double source_monomial_local[NALPHA] = {0};

		// if index in bounds
		if(idx < n){
			
			// calculate dx/h
			for(int d=0; d<DIMENSIONS; d++){
				dx[d] = (x[idx*DIMENSIONS+d] - centers[cidx[idx]*DIMENSIONS + d])/h;
			}
		    // calculate source monomial
			calc_monomial_gpu(idx, DIMENSIONS, pmax, nalpha, dx, source_monomial_local);
			dx2 = l2normsq_gpu(dx);
			ex = exp(-dx2)*q[idx];

			// source*exp(-dx^2/h^2),  2^{alpha}/alpha! is multiplied in the reduction step to save the computations
			for(int alpha=0; alpha<NALPHA; alpha++){
				c_alpha_dev[cidx[idx]*NALPHA*n + alpha*n + idx] = source_monomial_local[alpha]*ex;
			}
		}
		
	}
__global__ void reduce_k(int n, int nk, double *c_alpha_dev){
		extern __shared__ double tmp_c[];
		int tid = threadIdx.x;
		int blockdim = blockDim.x;
		int idx = tid + 2*blockdim*blockIdx.x;

		for(int k=0; k< nk; k++){
			for(int alpha=0; alpha<NALPHA; alpha++){
				if(idx < n){
					tmp_c[tid] = c_alpha_dev[k*NALPHA*n + alpha*n + idx];
				}
				else{
					tmp_c[tid] = 0.0;
				}
				if(idx + blockdim < n){
					tmp_c[tid + blockdim] = c_alpha_dev[k*NALPHA*n + alpha*n + idx + blockdim];
				}
				else{
					tmp_c[tid + blockdim] = 0.0;
				}
				__syncthreads();
				// if(fabs(tmp_c[tid]) > 1e-6){
				// 	printf("%f %f\n", tmp_c[tid], tmp_c[tid+blockdim]);
				// }
				// __syncthreads();
				for(int i=blockdim; i >= 1; i = i/2 ){
					if(tid < i){
						tmp_c[tid] += tmp_c[tid + i];
					}
					__syncthreads();
				}
				if(tid == 0){
					AtomicAdd_8(&c_alpha[k*NALPHA + alpha], tmp_c[tid]*constant[alpha]);
					//printf("%f %d\n", tmp_c[tid], k);
					//if(idx < 51000)
					//	printf("%d\n", idx);
					//c_alpha_dev[k*NALPHA*n + alpha*n + idx] = tmp_c[0];
					
				}
				__syncthreads();
			}
		}

	}


	__global__ void reduce_k2(int n, int nk, double *c_alpha_dev, int blocksize){
		extern __shared__ double tmp_c[];
		int tid = threadIdx.x;
		int blockdim = blockDim.x;
		int leftidx = 2*blocksize*tid;
		int rightidx = 2*blocksize*tid + blockdim*2*blocksize;
		printf("%d %d\n", leftidx, rightidx);
		for(int k=0; k < nk; k++){
			for(int alpha=0; alpha<NALPHA; alpha++){
				__syncthreads();
				if(leftidx < n){
					tmp_c[tid] = c_alpha_dev[k*NALPHA*n + alpha*n + leftidx];
				}
				else{
					tmp_c[tid] = 0.0;
				}
				if(rightidx < n){
					tmp_c[tid + blockdim] = c_alpha_dev[k*NALPHA*n + alpha*n + rightidx];
				}
				else{
					tmp_c[tid + blockdim] = 0.0;
				}
				__syncthreads();
				for(int i=blockdim; i >= 1; i = i/2 ){
					if(tid < i){
						tmp_c[tid] += tmp_c[tid + i];
					}
					__syncthreads();
				}
				if(tid == 0){
					c_alpha[k*NALPHA + alpha] = tmp_c[tid];
					//AtomicAdd_8(&c_alpha[k*NALPHA + alpha], tmp_c[tid]);
				}
				__syncthreads();
			}
		}

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
	}	void calc_centers_forgpu(int n, int nk, double *centers, unsigned int *cidx, int *delta_cidx, double *x, unsigned int *np, unsigned int *delta){
		unsigned int old_cidx;
		unsigned int new_cidx;

	*delta = 0;
		for(int i=0; i<n; i++){
			if(delta_cidx[i] != 0){

				old_cidx = cidx[i];
				new_cidx = old_cidx + delta_cidx[i];
				cidx[i] = new_cidx;
				assert(old_cidx + delta_cidx[i] >= 0);
				for(int d=0; d<DIMENSIONS; d++){
					if(np[old_cidx] - 1 > 0)
						centers[old_cidx*DIMENSIONS+d] = (np[old_cidx]*centers[old_cidx*DIMENSIONS+d] - x[i])/(np[old_cidx] - 1);
					else
						centers[old_cidx*DIMENSIONS+d] = 0;
					centers[new_cidx*DIMENSIONS+d] = (np[new_cidx]*centers[new_cidx*DIMENSIONS+d] + x[i])/(np[new_cidx] + 1);

				}
				np[old_cidx] -= 1;
				np[new_cidx] += 1;
				*delta += 1;

			}
		}

	}



//__global__ void clustering_copy_centers(int ncluster, double *clus_centers_dev, double *clus_centers_new_dev, unsigned int *clus_np_dev, unsigned int *clus_np_new_dev);
	// __global__ void clustering_swap(int ncluster, double *clus_centers_dev, double *clus_centers_new_dev, unsigned int *clus_np_dev, unsigned int *clus_np_new_dev);

	__global__ void clustering_copy_centers(int ncluster, double *centers_dev, double *centers_new_dev, unsigned int *np_dev, unsigned int *np_new_dev){
		int idx = blockDim.x*blockIdx.x + threadIdx.x;
		if(idx < ncluster*DIMENSIONS){
			centers_new_dev[idx] = centers_dev[idx];
		}
		if(idx < ncluster){
			np_new_dev[idx] = np_dev[idx];
		}
	}

	__global__ void clustering_swap(int ncluster, double *centers_dev, double *centers_new_dev, unsigned int *np_dev, unsigned int *np_new_dev){
		int idx = blockDim.x*blockIdx.x + threadIdx.x;
		if(idx < ncluster*DIMENSIONS){
			centers_dev[idx] = centers_new_dev[idx];
		}

		if(idx < ncluster){
			np_dev[idx] = np_new_dev[idx];
		}
	}