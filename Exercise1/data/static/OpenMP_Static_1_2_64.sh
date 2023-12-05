#!/bin/bash
#SBATCH --partition=EPYC 
#SBATCH --job-name=OS1264
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2 
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=OpenMP_Static_1_2_64.out


module load openMPI/4.1.5/gnu

MAPBY=socket
export OMP_PLACES=cores
export OMP_PROC_BIND=close


n=100

name=OpenMP_Static_1_2_64.csv
mode=static
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

for k in 10000 20000
do
    formatted_number=$(printf "%05d" "$k")
    filename="init_"$formatted_number".pgm" 
    for nthreads in 1 2 4 8 16 20 28 32 40 48 56 60 64
    do
        export OMP_NUM_THREADS=$nthreads
        for count in {1..5..1}
        do
            mpirun -n 2 --map-by $MAPBY  ./main.x -r -k $k -n $n -e $e -f $filename >> $res
        done
    done
done

