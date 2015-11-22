#
# Franklin and Hopper
#
CC = g++ -std=c++11
MPCC = g++ -std=c++11
OPENMP = -mp -fopenmp
CFLAGS = -O3
LIBS =

#
# Bassi  - Being retired soon/ Don't have access
#
#CC = cc -+ -qsuppress=1500-036
#MPCC = mpcc
#OPENMP = -qsmp=omp
#CFLAGS = -O3
#LIBS = -lm

#
# Jacquard
#
#CC = pathCC
#MPCC = mpicxx 
#OPENMP = -mp
#LIBS = -lm
#CFLAGS = -O3

#
# DaVinci - Being retired soon/ Don't have access
#
#CC = g++
#MPCC = g++
#OPENMP = -fopenmp
#LIBS = -lm
#CFLAGS = -O3
#MPILIBS = -lmpi++ -lmpi

TARGETS =  mpi

all:	$(TARGETS)

serial: serial.o common.o grid.o
	$(CC) -o $@ $(LIBS) serial.o common.o
pthreads: pthreads.o common.o grid.o
	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o
openmp: openmp.o common.o grid.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: 	mpi.o common.o grid.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o grid.o

openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp
mpi.o: mpi.cpp mpi.h common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp
grid.o: grid.cpp grid.h
	$(CC) -c $(CFLAGS) grid.cpp

clean:
	rm -f *.o $(TARGETS)
