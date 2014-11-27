CXX = nvcc
NCXX = nvcc
CFLAGS = -g -std=c++11
DEPS = utils.h ifgt.h io.h utils_gpu.h ifgt_gpu.h
OBJ = main.o utils.o ifgt.o io.o utils_gpu.o ifgt_gpu.o 
INC = -I/usr/local/cuda-6.5/include
%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(INC)

%.o: %.cu $(DEPS)
	$(NCXX) -ccbin=g++ -dc -o $@ $< $(INC) $(CFLAGS)


main: $(OBJ)
	$(NCXX) $^ $(CFLAGS) -dlink -o dlink.o $(INC)
	$(CXX) -o has.ty $^ dlink.o $(CFLAGS) -lcuda -lcudart -L/usr/local/cuda-6.5/lib64 $(INC)
	
clean:
	rm *.o
	rm has.ty