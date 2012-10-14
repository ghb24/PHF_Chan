SHELL = /bin/sh
CPP = g++
CSRCS = grid_3d_class.cpp main.cpp
COBJS = grid_3d_class.o main.o
PGM = test

$(PGM):	$(COBJS)
	$(CPP) -o $(PGM) $(COBJS)


grid_3d_class.o:	grid_3d_class.h grid_3d_class.cpp
	$(CPP) -c grid_3d_class.cpp

main.o:		grid_3d_class.h grid_3d_class.cpp main.cpp
	$(CPP) -c main.cpp

clean:
	rm *.o ; rm test;
