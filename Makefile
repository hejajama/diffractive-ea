#CXXFLAGS = `gsl-config --cflags` -g #-fopenmp #-pg -Wall
CXXFLAGS = `gsl-config --cflags` -O2 -I /nashome2/hejajama/amplitudelib/  
#CXXFLAGS = `gsl-config --cflags` -O2 -fopenmp
LDFLAGS = `gsl-config --libs` -fopenmp 

SOURCES = src/dipole.cpp src/gaus_lc.cpp src/vector.cpp src/nucleus.cpp \
	src/dipxs.cpp src/gdist.cpp \
	src/mersenne/mersenne_inline.cpp src/wave_function.cpp \
	src/dipxs_ipsat.cpp  \
	src/dipxs_ipnonsat.cpp src/dipxs_iim.cpp \
	src/gdist/gdist_dglap.cpp src/calculator.cpp \
	src/gauss_boost.cpp src/dipxs_bk.cpp \
	src/dipxs_ipsat2012.cpp src/virtual_photon.cpp \
	src/ipsat_mz/dipoleamplitude.cpp src/ipsat_mz/dglap_cpp/AlphaStrong.cpp \
	src/ipsat_mz/dglap_cpp/EvolutionLO.cpp

OBJECTS=$(SOURCES:.cpp=.o)

#OBJECTS  = src/main.o src/dipole.o
#src/main.o: src/main.cpp
#src/dipole.o: src/dipole.cpp

all: dipole 

dipole: $(OBJECTS) src/main.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) src/main.cpp src/libColorDipole/libraries/libColorDipole.a /nashome2/hejajama/amplitudelib/libamplitude.a -o dipole -lgfortran

libColorDipole: src/libColorDipole/src/*.f
	cd src/libColorDipole/ && make TEST_DIPOLE FC=gfortran
	#$(MAKE) -C src/libColorDipole/ TEST_DIPOLE FC=gfortran

tools: $(OBJECTS) tools/tws.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) tools/tws.cpp src/libColorDipole/libraries/libColorDipole.a /nashome2/hejajama/amplitudelib/libamplitude.a -o tws -lgfortran

generator: tools/create_amplitudelib_from_ipsat.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -I. tools/create_amplitudelib_from_ipsat.cpp src/libColorDipole/libraries/libColorDipole.a  /nashome2/hejajama/amplitudelib/libamplitude.a   -o dipxs_generator -lgfortran 

.cpp.o:
	g++ $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f dipole
