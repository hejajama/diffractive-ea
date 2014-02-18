#CXXFLAGS = `gsl-config --cflags` -g #-fopenmp #-pg -Wall
CXXFLAGS = `gsl-config --cflags` -O2 
#CXXFLAGS = `gsl-config --cflags` -O2 -fopenmp
LDFLAGS = `gsl-config --libs` 

SOURCES = src/dipole.cpp src/gaus_lc.cpp src/vector.cpp src/nucleus.cpp \
	src/dipxs.cpp src/gdist.cpp \
	src/mersenne/mersenne_inline.cpp src/wave_function.cpp \
	src/dipxs_ipsat.cpp  \
	src/dipxs_ipnonsat.cpp src/dipxs_iim.cpp \
	src/gdist/gdist_dglap.cpp src/calculator.cpp \
	src/gauss_boost.cpp

OBJECTS=$(SOURCES:.cpp=.o)

#OBJECTS  = src/main.o src/dipole.o
#src/main.o: src/main.cpp
#src/dipole.o: src/dipole.cpp

all: dipole 

dipole: $(OBJECTS) src/main.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) src/main.cpp -o dipole 

generator:
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) create_amplitudelib_from_ipsat.cpp -o dipxs_generator 

.cpp.o:
	g++ $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f dipole
