#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=provampi
#SBATCH --nodes=1

#SBATCH --ntasks-per-node=2

#SBATCH --cpus-per-task=3
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=prova_mpi_global.out

##mpicc -fopenmp prova_mpi_global.c -o prova_mpi_global
##module load openMPI/4.1.5/gnu/12.2.1 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpirun -np 2 ./prova_mpi_global
