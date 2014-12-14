CXX = nvcc
NCXX = nvcc
CFLAGS = -std=c++11 -ccbin=g++ --maxrregcount=32 -lineinfo -lcurand 
DEPS = utils_host.h io.h utils_dev.h ifgt_dev.h clustering_dev.h ifgt_host.h ifgt_dev_utils.h clustering_host.h exact_gt.h
OBJ = main.o utils_host.o io.o utils_dev.o ifgt_dev.o clustering_dev.o ifgt_host.o ifgt_dev_utils.o clustering_host.o exact_gt.o
INC = -I/usr/local/cuda-6.5/include

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(INC) $(CFLAGS) 

%.o: %.cu $(DEPS)
	$(NCXX) -dc -o $@ $< $(INC) $(CFLAGS) 


main: $(OBJ)
	$(NCXX) $^ $(CFLAGS) -dlink -o dlink.o $(INC) 
	$(CXX) -o has.ty $^ dlink.o -lcuda -lcudart -L/usr/local/cuda-6.5/lib64 $(INC) $(CFLAGS)
		
clean:
	rm *.o
	rm has.ty
