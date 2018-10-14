main:
	@echo " Compile mpi ..."
	mpicc -o mpi ./mpi.c -lm
openMp:
	@echo " Compile openmp ..."
	mpicc -openmp -o openmp ./openmp.c -lm
runMpi:
	mpiexec -n 25 -f machines ./mpi 1920 2520  waterfall_1920_2520.raw 25 COLOR
runOpenMp:
	mpiexec -n 25 -f machines ./openmp 1920 2520 waterfall_1920_2520.raw 25 COLOR
clean:
	rm -f ./mpi ./openmp
