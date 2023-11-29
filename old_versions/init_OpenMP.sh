#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=OMP_Init
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
##SBATCH --mem=200gb
#SBATCH --mem=2gb
#SBATCH --time=00:30:00 
##SBATCH --exclusive
#SBATCH --output=OMP_Init.out

module load openMPI/4.1.5/gnu


# Define MPI binding and OMP affinity

export OMP_PLACES=cores
export OMP_PROC_BIND=close

init=data/init_OpenMP.csv



echo "size,threads,time" > $init

for ksize in 200 300
do
    formatted_number=$(printf "%05d" "$ksize")
    filename="init_"$formatted_number".pgm" 
    echo $filename
    for n_threads in 1 2 4 
    do
        export OMP_NUM_THREADS=$n_threads
        
        mpirun -np 1  ./main.x -i -k $ksize -f $filename > init.txt
        pwd
        time_value=$(grep -o 't_init: [0-9.]*' init.txt | awk '{print $3}')
        echo time_value: $time_value
        echo "$ksize,$n_threads,$time_value" >> $init

        rm init.txt
       
    done
done

##rm output_initialization_openMP.txt # Remove useless temporary file


##--map-by $MAPBY --bind-to $BINDTO