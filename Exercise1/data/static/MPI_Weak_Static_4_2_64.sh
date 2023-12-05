#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=MWA4264
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --mem=200gb 
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=MPI_Weak_Static_4_2_64.out

module load openMPI/4.1.5/gnu

MAPBY=node
BINDTO=socket
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=64

n=100
name=MPI_Weak_Static_4_2_64
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

ki=5000
for nprocs in 1 2 3 4 5 6 7 8 
do
    
    k=$(echo "scale=0; sqrt($nprocs * $ki^2)" | bc)
    formatted_number=$(printf "%05d" "$k")
    filename="init_"$formatted_number".pgm" 
   
    for rep in {1..5..1}
    do
            mpirun -n $nprocs --map-by $MAPBY --bind-to $BINDTO ./main.x -r -k $k -n $n -e $e -f $filename >> $res
    done
done




   






        