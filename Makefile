

OFILES = \
	 base/Geometry.o\
	 base/SPHSolver.o\
	 base/Window.o\




ROOTDIR = .
LIB = $(ROOTDIR)/lib/libstarter.a 


GLLDFLAGS     = -lglut -lGL -lm -lGLU


CXX = g++ -Wall -g -O2 -fPIC $(DEFINES) -fopenmp -std=c++11



INCLUDES =  -I ./include/ -I /usr/local/include/ -I/usr/include/ 




.C.o:
	$(CXX) -c $(INCLUDES) $< -o $@

base: $(OFILES)
	ar rv $(LIB) $?
	$(CXX) base/main.C $(INCLUDES)  -L./lib -lstarter $(GLLDFLAGS)  -o bin/sph

clean:
	rm -rf bin/sph  *.o base/*.o base/*~ include/*~  $(LIB)  




