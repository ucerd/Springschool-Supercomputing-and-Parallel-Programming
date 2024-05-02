# springschool24
[https://www.pakistansupercomputing.com/school.php]

Supercomputing and Parallel Programming Spring School

The Role and Importance of Supercomputing Centers in Third-Generation Entrepreneurial Institutes Understanding Supercomputing and HPC Architectures:
Learning Objective: understand software stack and hardware configuration of HPC cluster 
1. HPC Server (GP-GPU based System)
2. HPC Cluster (Computer Nodes, Parallel Programming Models and Libraries, Network Configuration, File System)
3. Supercomputing System Software Stack and Parallel Programming Libraries
4. Open Hardware and Software-Based Supercomputing Infrastructure


**C Programming for HPC Applications**
Learning Objective: understanding complicated codebases
1. Memory Addressing in Parallel Architectures
2. Dynamic Function Dispatching in Parallel Programming
3. Tracing impactful projects like Linux Kernel or OpenMPI


**Hands-on Accessing a Supercomputer**
Learning Objective: get comfortable with a supercomputing environment
1. ssh, scp
2. Slurm Commands: sinfo, sbatch, srun, squeue
3. Development in Linux Command Line: vim, tmux, gcc/clang
4. Linux Environment for Developers: PATH, LD_LIBRAY_PATH, L_PRELOAD


**Art of Parallel Programming: Think Parallel**
Learning Objective: understand code architecture, data dependencies, control and data flows
1. Understanding Parallelism
2. Parallel Programming Models
3. Designing Parallel Algorithms
4. Optimizing Parallel Performance


**Multicore Programming**
Learning Objective: understand why we program multithreaded codes and do it hands-on
1. Multicore Architecture
   a) What does a modern computer look like? Multicore, Vector architecture, (b) Processes and Threads, (c) How processes and threads are spawned in C
2. Shared Memory Programming Model: OpenMP
   b) OpenMP Model - shared memory (b) OpenMP Standard and features (c) OpenMP Implementations: Compilation and Runtime
3. OpenMP Programming
   c) parallel for b) SIMD c) task
OpenMP Tasking

Code: https://github.com/OpenMP/Examples

https://www.openmp.org/resources/refguides/

https://www.openmp.org/wp-content/uploads/openmp-examples-5.2.2-final.pdf


**The Current Mess (and Power) of Heterogeneity in HPC**
Multinode Programming: Learning Objective: understand why we program multinode codes and do it hands-on
1. Multinode Architecture
   a) Components of a supercomputer: login and compute nodes, offload devices, memory, storage and interconnect
   b) Interconnects: ethernet, infiniband, nvlink
   c) Network Technologies: Ethernet and Infiniband
2. Difference between ethernet and infiniband
   a ) Why is infiniband needed?
   b) How infiniband (RDMA) programming is different from sockets
   c) A glimpse of RDMA
   d) RDMA Architecture
3. RDMA Programming and why you should never do it

**Multinode Programming- Continue**
Learning Objective: understand why we program multinode codes and do it hands-on
1. MPI Programming
   a) What does MPI encapsulate?
   b) How does MPI do it? An OpenMPI perspective
   c) MPI Syntax: Send/Recv, ISend/IRecv, Wait
   d) MPI Collectives: allreduce, allgather, allscatter, all-to-all
   e) Further MPI Capabilities and Applications
   f) GPU-enabled MPI
   g) DPU-enabled MPI
   h) Advanced Data Types as Cartesian

Info: https://mpitutorial.com/tutorials/mpi-hello-world/


**Offload Programming**
Learning Objective: understand why we program offload codes and do it hands-on
1. Intro to Offloading Devices: GPUs, DPUs and FPGAs
   a) Offload Programming Approaches: Libraries, Directives & Domain Specific Languages
2. Directive-based Approach: OpenMP and OpenACC Offloading
   a) OpenMP Offload Architecture
   b) OpenMP Offload Features
   c) OpenMP GPU Offloading
   d) OpenMP DPU Offloading


**Distributed Artificial Intelligence** 
1. Understanding AI Models
2. Types of AI Model Programming
   a) Data Level Parallel and Model Level Parallel
   b) Single GPU based Execution
3. Multi-GPU based Execution
4. Multi-Node based Execution


**Videos Lectures @**

https://www.youtube.com/playlist?list=PLxtSKP68m3fo2qQmkNvRu4KPLbRb-UHla

