# define default variables
procs?=3
ruut?=1


comp:
	mpiCC -o project_1_2 project_1_2.cpp utils.cpp

run:
	mpi_run -np $(procs) ./project_1_2 $(ruut)

debug:
	mpiCC -g project_1_2.cpp

gdb:
	mpiexec -d -np $(procs) ./project_1_2 $(ruut)
	# mpiexec -n $(procs) gdb ./a.out

clean:
	rm project_1_2 slurm*