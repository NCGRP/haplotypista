# serial compilation
CXX = g++
FLAGS = -O3 -Wall

# for compilation on Linux and Mac desktops:
# 1. use mpic++ compiler, and -fopenmp directive
#CXX = mpic++
#OMP = -fopenmp

# for compilation on CSU ISTEC Cray 
# 1. make clean
# 2. use CC compiler, and -fopenmp directive
#CXX = CC
#OMP = -fopenmp

# for compilation on TACC Lonestar
# 1. make clean
# 2. use mpicxx compiler, and -openmp directive
#CXX = mpicxx
#OMP = -openmp

#FLAGS = -O3 -Wall

all: h.o
	${CXX} ${FLAGS} h.o -o haplotypista

h.o: h.cpp h.hpp
	${CXX} ${FLAGS} -c h.cpp 


clean:
	rm -rf *.o
