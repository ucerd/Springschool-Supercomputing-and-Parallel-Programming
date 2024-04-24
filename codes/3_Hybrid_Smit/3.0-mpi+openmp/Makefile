include ../Makefile.common

TASK_PATH = $(PWD)
DEBUG =
CFLAGS_EXTRA = -fopenmp $(DEBUG)
OMP_NUM_THREADS = 2
CC = mpicc
BLOCK = 2
CHUNK = 2

mpi: sw
	mpirun -np 2 ./sw $(DATAPATH)/a_500k $(DATAPATH)/b_500k \
	                  $(DATAPATH)/data.score $(PENALTY) $(SIZE) \
	                  $(BLOCK) $(BLOCK) $(CHUNK)
