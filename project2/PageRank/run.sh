#!/bin/bash

if [ $1 == 0 ]
then
    time ./PageRank $2 $3 > PageRank.$2.output
    tail -1 PageRank.$2.output > output.txt
else
    ./PageRank $2 $3
fi