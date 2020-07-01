#!/bin/bash

if [ $2 == 0 ]
then
    mpirun -np $1 ./MPI-WordsCount -f ./big_file/big_100.txt
else
    mpirun -np $1 ./MPI-WordsCount -d ./small_file/
fi