#!/bin/bash
echo "test 1000 times by using 10 threads"
./mcopm 100 10
echo "test 1000 times by using 5  threads"
./mcopm 200 5
echo "test 1000 times by using 1 threads"
./mcopm 1000 1
echo "test 100000000 times by using 10 threads"
./mcopm 10000000 10
echo "test 100000000 times by using 5  threads"
./mcopm 20000000 5
echo "test 100000000 times by using 1  threads"
./mcopm 100000000 1