#CXXFLAGS = `gsl-config --cflags` -g -pg 
CXXFLAGS = `gsl-config --cflags` -O2 
#CXXFLAGS = `gsl-config --cflags` -O2 -fopenmp
LDFLAGS = `gsl-config --libs` 

SOURCES = main.cpp gdist_dglap.cpp
OBJECTS=$(SOURCES:.cpp=.o)

#OBJECTS  = src/main.o src/dipole.o
#src/main.o: src/main.cpp
#src/dipole.o: src/dipole.cpp

all: gdist 

gdist: $(OBJECTS)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) -o gdist

.cpp.o:
	g++ $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS)
	rm -f dipole
