#!/bin/bash
#SBATCH --partition=EPYC 
#SBATCH --job-name=OMP_Static
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2 
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=OpenMP_scalability2.out


module load openMPI/4.1.5/gnu

MAPBY=node
BINDTO=socket
export OMP_PLACES=cores
export OMP_PROC_BIND=close

name=$(basename $0 .sh)

echo $name

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
            for nthreads in 28 32 38 42 48 52 58 62 64
            do
                export OMP_NUM_THREADS=$nthreads
                for count in {1..5..1}
                do
                    mpirun -n 4 --map-by $MAPBY --bind-to $BINDTO ./main.x -r -k $k -n $n -e $e -f $filename >> $res
                done
            done
        done
    elif [ $e -eq 0 ]; then
        k=10000
        formatted_number=$(printf "%05d" "$k")
        filename="init_"$formatted_number".pgm" 
        for nthreads in 28 32 38 42 48 52 58 62 64
        do
            export OMP_NUM_THREADS=$nthreads
            for count in {1..5..1}
            do
                mpirun -n 4 --map-by $MAPBY --bind-to $BINDTO ./main.x -r -k $k -n $n -e $e -f $filename >> $res
            done
        done
    fi
done