#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=main3
#SBATCH --nodes=1

#SBATCH --ntasks-per-node=2

#SBATCH --cpus-per-task=3
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=main3.out

##mpicc -fopenmp prova_mpi_global.c -o prova_mpi_global
##module load openMPI/4.1.5/gnu/12.2.1 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np 2 ./main3.x -f init_00100.pgm -r -e 1 -n 10 -s 1 -k 100 
