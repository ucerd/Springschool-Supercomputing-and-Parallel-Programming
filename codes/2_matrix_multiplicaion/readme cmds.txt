export PATH="/usr/local/cuda/bin:$PATH"
gcc -o matmul_x86 ./matmul.c
gcc -fopenmp -o multi_core ./matmul_omp.c 
export OMP_NUM_THREADS=4
export OMP_NUM_THREADS=1
nvcc -arch=sm_86 -o cuda_matmul matmul.cu 


#master MPI
mpicc matmul240.c -o matmul240
mpiexec --oversubscribe -np 240 ~/work/code/mpi/matmul/matmul240


#SCP
scp ./matmul240 node01:/home/caid_main/work/code/mpi/matmul/
scp ./matmul240 node02:/home/caid_main/work/code/mpi/matmul/


#Multi-Node
mpiexec -np 240 --oversubscribe --host master,node01,node02 ~/work/code/mpi/matmul/matmul240
