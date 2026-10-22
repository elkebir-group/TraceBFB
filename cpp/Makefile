CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDFLAGS ?= -pthread

OBJECTS = nbfb_cli.o nbfb.o nbfb_blocks.o nbfb_simple.o nbfb_state.o

all: nbfb

nbfb: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

nbfb_cli.o: nbfb_cli.cpp nbfb.hpp nbfb_blocks.hpp nbfb_simple.hpp nbfb_state.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_cli.cpp

nbfb.o: nbfb.cpp nbfb.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb.cpp

nbfb_blocks.o: nbfb_blocks.cpp nbfb_blocks.hpp nbfb.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_blocks.cpp

nbfb_simple.o: nbfb_simple.cpp nbfb_simple.hpp nbfb.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_simple.cpp

nbfb_state.o: nbfb_state.cpp nbfb_state.hpp nbfb.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_state.cpp

clean:
	rm -f nbfb $(OBJECTS)

.PHONY: all clean
