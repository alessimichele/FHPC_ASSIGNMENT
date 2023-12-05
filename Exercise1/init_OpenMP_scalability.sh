#!/bin/bash
#SBATCH --partition=EPYC 
#SBATCH --job-name=init_OpenMP_scalability
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=init_OpenMP_sacalability.out

module load openMPI/4.1.5/gnu
MAPBY=node
BINDTO=socket
export OMP_PLACES=cores
export OMP_PROC_BIND=close

name=$(basename $0 .sh)
name=${name#init_}
echo $name

mode=init

res=../data/$mode/$name.csv

echo "mode,size,nthreads,k,time" > $res

for k in 10000 20000
do
    formatted_number=$(printf "%05d" "$k")
    filename="init_"$formatted_number".pgm" 
    for nthreads in 1 2 4 8 16 22 28 32 38 42 48 52 58 64
    do
        export OMP_NUM_THREADS=$nthreads
        for rep in {1..5..1}
        do
            mpirun -n 4 --map-by $MAPBY --bind-to $BINDTO ../main.x -i -k $k -f $filename >> $res
        done
    done
done
