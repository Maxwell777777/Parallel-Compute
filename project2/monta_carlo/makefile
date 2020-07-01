EXECS=mcopm
MPICC?=mpicc

all: ${EXECS}

mcopm: mcopm.o
	${MPICC} mcopm.o -lm -o mcopm
    
mcopm.o: mcopm.c
	${MPICC} -c mcopm.c
  
clean:
	rm -f ${EXECS}
