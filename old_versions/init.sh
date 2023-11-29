#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=OMP_Init
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --mem=500Mb
#SBATCH --time=00:05:00 

#SBATCH --output=OMP_Init.out

module load openMPI/4.1.5/gnu


# Define MPI binding and OMP affinity
MAPBY=node
BINDTO=socket

export OMP_PLACES=cores
export OMP_PROC_BIND=close

init=data/init_OpenMP.csv

echo "size,threads,time" > $init

ksize=20
n_threads=8

formatted_number=$(printf "%05d" "$ksize")

filename="init_"$formatted_number".pgm" 

export OMP_NUM_THREADS=$n_threads

mpirun -n 2 ./main.x -i -k $ksize -f $filename 



       
    


##rm output_initialization_openMP.txt # Remove useless temporary file

## time_value=$(grep -o 't_init: [0-9.]*' init.txt | awk '{print $3}')
##--map-by $MAPBY --bind-to $BINDTO