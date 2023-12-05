#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=MPI_Strong
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16 
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=MPI_strong1.out

module load openMPI/4.1.5/gnu

name=$(basename $0 .sh)

MAPBY=numa
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=16

n=100

for mode in static wave ordered
do
    res=data/$mode/$name.csv
    if [ "$(ls -A data/$mode/)" ]; then
        truncate -s 0 $res
    else
        echo "mode,size,nthreads,k,time" > $res
    fi

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

    if [ $e -eq 1 ] || [ $e -eq 2 ]; then
        for k in 10000 20000
        do
            formatted_number=$(printf "%05d" "$k")
            filename="init_"$formatted_number".pgm" 
            for nprocs in 1 2 4 6 8 
            do
                for rep in {1..5..1}
                do
                    mpirun -n $nprocs --map-by $MAPBY ./main.x -r -k $k -n $n -e $e -f $filename >> $res
                done
            done
        done
    elif [ $e -eq 0 ]; then
        k=10000
        formatted_number=$(printf "%05d" "$k")
        filename="init_"$formatted_number".pgm" 
        for nprocs in 1 2 4 6 8 
        do
            for rep in {1..5..1}
            do
                mpirun -n $nprocs --map-by $MAPBY ./main.x -r -k $k -n $n -e $e -f $filename >> $res
            done
        done
    fi
done





        