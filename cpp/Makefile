CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDFLAGS ?= -pthread

all: nbfb

nbfb: nbfb_cli.o nbfb_fast.o nbfb_state.o nbfb_blocks.o
	$(CXX) $(CXXFLAGS) -o $@ nbfb_cli.o nbfb_fast.o nbfb_state.o nbfb_blocks.o $(LDFLAGS)

nbfb_cli.o: nbfb_cli.cpp nbfb_fast.hpp nbfb_state.hpp nbfb_blocks.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_cli.cpp

nbfb_fast.o: nbfb_fast.cpp nbfb_fast.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_fast.cpp

nbfb_state.o: nbfb_state.cpp nbfb_state.hpp nbfb_fast.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_state.cpp

nbfb_blocks.o: nbfb_blocks.cpp nbfb_blocks.hpp nbfb_fast.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_blocks.cpp

clean:
	rm -f nbfb nbfb_cli.o nbfb_fast.o nbfb_state.o nbfb_blocks.o

.PHONY: all clean
