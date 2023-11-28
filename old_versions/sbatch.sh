#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=prova
#SBATCH --nodes=1

#SBATCH --ntasks-per-node=2

#SBATCH --cpus-per-task=4

#SBATCH --time=00:10:00 

#SBATCH --output=prova.out

module load openMPI/4.1.5/gnu/12.2.1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np 2 ./main.x -f init_00010.pgm -r -e 2 -n 10 -s 2 -k 10

##mpirun -np 2 ./main.x -f init_00010.pgm -i -k 10


##mpicc -fopenmp prova_mpi_global.c -o prova_mpi_global
##module load openMPI/4.1.5/gnu/12.2.1 
## mpirun -np 2 ./main.x -f init_00010.pgm -r -e 0 -n 10 -s 2 -k 10

