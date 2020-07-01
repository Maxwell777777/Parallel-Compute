Implement the matrix multiplication function, Use a 4 * 4 kernel to do the pooling and convolution operation for the matrix.
1. Compile the source code in ./covolution and ./supporting, the makefile is provided. Move project_1_2 to ./.
2. Compile the main function: mpicc -O -o main main.c -lm.
3. ./supporting/seqmm mat1 mat2: sequential matmul;
    ./supporting/prtmat mat x: print the first x rows and columns of mat;
    ./supporting/genmat x1 x2 mat: generate a matrix with x1 rows and x2 columns;
4. Run run.sh to do the matmul. There are 2 parameters, 1) the dimension of the matrix; 2) the number of threads.
5. Run conv.sh to do pooling and convolution. The parameter is matrix filename.
6. Run run_cluster.sh to do the MPI_CLUSTER matmul. There are 2 parameters, 1) the dimension of the matrix; 2) the number of threads.