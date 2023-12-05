#!/bin/bash
#SBATCH --partition=EPYC 
#SBATCH --job-name=OW1164
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64 ## 128 if MAPBY=node
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=OpenMP_Wave.out


module load openMPI/4.1.5/gnu

MAPBY=socket # node
export OMP_PLACES=sockets # cores threads
export OMP_PROC_BIND=false # close spread


echo $name

n=100

## OpenMP_Wave_proc_places_bind

name=OpenMP_Wave_s_sf
mode=wave
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

for k in 5000 
do
    formatted_number=$(printf "%05d" "$k")
    filename="init_"$formatted_number".pgm" 
    for nthreads in 1 8 16 24 32 40 48 56 64 2 3 4 6 12 20 28 36 4 50 58
    do
        export OMP_NUM_THREADS=$nthreads
        for count in {1..5..1}
        do
            mpirun -n 1 --map-by $MAPBY ./main.x -r -k $k -n $n -e $e -f $filename >> $res
        done
    done
done

