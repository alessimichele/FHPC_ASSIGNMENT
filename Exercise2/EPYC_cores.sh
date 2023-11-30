#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=bliscore
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=128
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=cores_epyc_blis.out

module load mkl
module load openBLAS/0.3.23-omp
export LD_LIBRARY_PATH=/u/dssc/mcarol00/blis/myblis/lib:$LD_LIBRARY_PATH

export OMP_PLACES=cores
export OMP_PROC_BIND=false #automatic binding


name=cores_scalability_study_EPYC_blis
res=EPYC/data/$name.csv
echo "library,precision,n_threads,iteration,time,GFLOPS" > $res

srun -n1 make cpu # Now I have all the needed executables.

k=10000

for library in 'blis' #'mkl' 'blis' 'oblas'
do
    for precision in 'f' 'd' # Single and double precision
    do
        for n_threads in 1 2 {4..128..4}
        do 
            for j in 1 2 3 4 5
            do
                export OMP_NUM_THREADS=$n_threads
                export BLIS_NUM_THREADS=$n_threads
                echo -n "$library,$precision,$n_threads,$j," >> $res
                srun -n1 --cpus-per-task=$n_threads ./"$precision"gemm_"$library".x $k $k $k >> $res

            done
        done
    done
done