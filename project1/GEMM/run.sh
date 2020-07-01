#!/bin/bash


./supporting/genmat $1 $1 mat$1.1
./supporting/genmat $1 $1 mat$1.2
echo "mat1:"
if [ $1 -gt 10 ]
then
    echo "the first 10 rows and cols"
    ./supporting/prtmat mat$1.1 10
else    
    ./supporting/prtmat mat$1.1
fi
echo "mat2:"
if [ $1 -gt 10 ]
then
    echo "the first 10 rows and cols"
    ./supporting/prtmat mat$1.2 10
else
    ./supporting/prtmat mat$1.2
fi
mpirun -np $2 ./main mat$1.1 mat$1.2 output_mat$1
echo "result:"
if [ $1 -gt 10 ]
then
    echo "the first 10 rows and cols"
    ./supporting/prtmat output_mat$1 10
else
    ./supporting/prtmat output_mat$1
fi
./supporting/seqmm mat$1.1 mat$1.2 output_mat_se$1
echo "sequential result:"
if [ $1 -gt 10 ]
then
    echo "the first 10 rows and cols"
    ./supporting/prtmat output_mat_se$1 10
else
    ./supporting/prtmat output_mat_se$1
fi