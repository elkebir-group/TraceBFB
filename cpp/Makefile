CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
LDFLAGS ?= -pthread

all: nbfb

nbfb: nbfb_cli.o nbfb_fast.o
	$(CXX) $(CXXFLAGS) -o $@ nbfb_cli.o nbfb_fast.o $(LDFLAGS)

nbfb_cli.o: nbfb_cli.cpp nbfb_fast.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_cli.cpp

nbfb_fast.o: nbfb_fast.cpp nbfb_fast.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ nbfb_fast.cpp

clean:
	rm -f nbfb nbfb_cli.o nbfb_fast.o

.PHONY: all clean
