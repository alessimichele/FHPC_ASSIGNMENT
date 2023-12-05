#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=MSO
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16 
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=MPI_Strong_Ordered.out

module load openMPI/4.1.5/gnu


MAPBY=numa
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16

n=100
name=MPI_Strong_Ordered
mode=ordered
res=data/$mode/$name.csv


if [ "$mode" == "static" ]
then
    e=1
elif [ "$mode" == "ordered" ]
then
    e=0
elif [ "$mode" == "wave" ]
then
    e=2
else
    echo "Error: mode not recognized"
    exit 1
fi

for k in 200 300 3000 5000 10000 
do
    formatted_number=$(printf "%05d" "$k")
    filename="init_"$formatted_number".pgm" 
    for nprocs in 1 2 4 6 8 10 12 14 16
    do
        for rep in {1..5..1}
        do
            mpirun -n $nprocs --map-by $MAPBY ./main.x -r -k $k -n $n -e $e -f $filename >> $res
        done
    done
done
   
   






        