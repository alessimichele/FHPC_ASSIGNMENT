#!/bin/bash
#SBATCH --partition=EPYC 
#SBATCH --job-name=IS
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=init_MPI_strong.out

module load openMPI/4.1.5/gnu

MAPBY=numa
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16



name=init_MPI_Strong
mode=init
res=data/$mode/$name.csv

for k in 10000 20000
do
    formatted_number=$(printf "%05d" "$k")
    filename="init_"$formatted_number".pgm" 
    for nprocs in 1 2 4 8 10 12 14 16 
    do
        for rep in {1..5..1}
        do
            mpirun -n $nprocs --map-by $MAPBY ./main.x -i -k $k -f $filename >> $res
        done
    done
done
