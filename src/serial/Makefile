CXX = g++
CFLAGS = -Wall -g -std=c++11
DEPS = utils.h ifgt.h io.h
OBJ = main.o utils.o ifgt.o io.o

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) 

main: $(OBJ)
	$(CXX) -o has.ty $^ $(CFLAGS)
	rm *.o

clean:
	rm has.ty