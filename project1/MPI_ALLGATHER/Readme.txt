Use MPI_SEND and MPI_RECV to implement the MPI_ALLGATHER function.
1) Compile the source code: mpicc -O -o ALLGATHER ALLGATHER.c.
2) Run run.sh. There are two parameters: 1) the number of threads; 2) the number of random numbers.
3) MPI Cluster: mpirun -n 8 -f /home/zq/MJT/MPI_ALLGATHER/config /home/zq/MJT/MPI_ALLGATHER/ALLGATHER 100000