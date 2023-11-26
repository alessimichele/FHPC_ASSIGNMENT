#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=provampi
#SBATCH --nodes=1
##i commenti sono scritti con due cancelletti, ntasks-per-node è il numero di processi per nodo ->mpi
#SBATCH --ntasks-per-node=2
##cpus-per-task è il numero di core (thread) per processo ->openmp
#SBATCH --cpus-per-task=3
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=prova_mpi_global.out

##mpicc -fopenmp prova_mpi_global.c -o prova_mpi_global

export OMP_NUM_THREADS=8
mpirun -np 2 ./prova_mpi_global
