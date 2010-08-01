CXXFLAGS = `gsl-config --cflags` -g
LDFLAGS = `gsl-config --libs` 

SOURCES = src/main.cpp src/dipole.cpp src/vm_photon.cpp src/vector.cpp src/nucleus.cpp \
	src/dipxs_ipnonsat.cpp src/dipxs.cpp src/gdist.cpp src/mersenne/mersenne_inline.cpp
OBJECTS=$(SOURCES:.cpp=.o)

#OBJECTS  = src/main.o src/dipole.o
#src/main.o: src/main.cpp
#src/dipole.o: src/dipole.cpp

all: dipole 

dipole: $(OBJECTS)
	g++ $(LDFLAGS) $(OBJECTS) -o dipole 

.cpp.o:
	g++ $(CFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f dipole
