Implement the wordcount algorithm with MPI.
1) mkdir build
2) cd build
3) cmake ..
4) make
5) mv MPI-WordsCount ..
6) cd ..
7) Run run.sh. There are two parameters: 1) the number of threads; 2) mode, if mode == 0, count big file; else count small file list.
8) MPI_CLUSTER: 
mpirun -n 8 -f /home/zq/MJT/WordCount/config /home/zq/MJT/WordCount/MPI-WordsCount -f /home/zq/MJT/WordCount/big_file/big_100.txt
mpirun -n 8 -f /home/zq/MJT/WordCount/config /home/zq/MJT/WordCount/MPI-WordsCount -d /home/zq/MJT/WordCount/small_file/
