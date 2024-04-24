include ../Makefile.common

TASK_PATH = $(PWD)
DEBUG =
CFLAGS_EXTRA = -fopenmp $(DEBUG)
OMP_NUM_THREADS = 1
CC = mpicc
BLOCK = 2

mpi: sw
	mpirun -np 1 ./sw $(DATAPATH)/a_500k $(DATAPATH)/b_500k \
	                  $(DATAPATH)/data.score $(PENALTY) $(SIZE) $(BLOCK) $(BLOCK)
