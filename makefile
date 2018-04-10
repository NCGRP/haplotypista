# serial compilation
CXX = g++
FLAGS = -std=c++11 -O3 -Wall

all: h.o
	${CXX} ${FLAGS} h.o -o haplotypista

h.o: h.cpp h.hpp
	${CXX} ${FLAGS} -c h.cpp 


clean:
	rm -rf *.o
