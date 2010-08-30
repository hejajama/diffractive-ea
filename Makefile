CXXFLAGS = `gsl-config --cflags` -g #-fopenmp #-pg -Wall
#CXXFLAGS = `gsl-config --cflags` -O2 
#CXXFLAGS = `gsl-config --cflags` -O2 -fopenmp -g
LDFLAGS = `gsl-config --libs` 

SOURCES = src/main.cpp src/dipole.cpp src/vm_photon.cpp src/vector.cpp src/nucleus.cpp \
	src/dipxs.cpp src/gdist_toy.cpp \
	src/mersenne/mersenne_inline.cpp src/wave_function.cpp \
	src/dipxs_ipsat.cpp src/dipxs_fourier.cpp \
	src/dipxs_ipnonsat.cpp src/dipxs_iim.cpp \
	src/gdist/gdist_dglap.cpp src/calculator.cpp

OBJECTS=$(SOURCES:.cpp=.o)

#OBJECTS  = src/main.o src/dipole.o
#src/main.o: src/main.cpp
#src/dipole.o: src/dipole.cpp

all: dipole 

dipole: $(OBJECTS)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o dipole 

.cpp.o:
	g++ $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f dipole
