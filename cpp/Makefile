CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDFLAGS ?= -pthread

OBJECTS = nbfb_cli.o nbfb.o nbfb_simple.o

all: nbfb

nbfb: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

nbfb_cli.o: nbfb_cli.cpp nbfb.hpp nbfb_simple.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_cli.cpp

nbfb.o: nbfb.cpp nbfb.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb.cpp

nbfb_simple.o: nbfb_simple.cpp nbfb_simple.hpp nbfb.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_simple.cpp

clean:
	rm -f nbfb $(OBJECTS)

.PHONY: all clean
