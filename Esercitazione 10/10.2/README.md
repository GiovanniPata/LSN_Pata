# ISTRUZIONI

- _COMPILAZIONE:_

		module load mpi/mpich-x86_64

		make

		mpicxx main.cpp
		


- _ESECUZIONE:_

		./mpiexec -np 8 a.out

- _OSSERVAZIONI:_ Necessario avere una installazione di MPI. Programma settato con nodi comunicanti. Per non farli comunicare è sufficiente alzare il valore di G_migr finché non diventa più alto del numero G di generazioni