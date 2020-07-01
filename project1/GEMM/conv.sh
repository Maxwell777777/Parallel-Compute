#!/bin/bash

if [ $2 == 0 ]
then
    echo "pooling"
    mpirun -np 4 ./project_1_2 $1 test_kernel_4_4_avg 0
else
    echo "conv"
    mpirun -np 4 ./project_1_2 $1 kernel 0
fi