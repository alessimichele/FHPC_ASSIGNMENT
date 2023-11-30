#!/bin/bash 
#SBATCH --partition=EPYC 
#SBATCH --job-name=oblasmkl
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=128
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=oblasmkl.out

module load mkl
module load openBLAS/0.3.23-omp
export LD_LIBRARY_PATH=/u/dssc/mcarol00/blis/myblis/lib:$LD_LIBRARY_PATH

export OMP_PLACES=cores
export OMP_NUM_THREADS=64
export BLIS_NUM_THREADS=64

name=complete_scalability_study_EPYC_oblas_mkl
res=EPYC/data/$name.csv
echo "library,precision,k,policy,iteration,time,GFLOPS" > $res

srun -n1 make cpu # Now I have all the needed executables.

for library in 'oblas' 'mkl' #'mkl' 'blis' 'oblas'
do
    for precision in 'f' 'd' # Single and double precision
    do
        for k in {2000..20000..1000} # Matrix size
        do
            for policy in 'false' 'close' 'spread'
            do
                for j in 1 2 3 4 5
                do
                    export OMP_PROC_BIND=$policy
                    echo -n "$library,$precision,$k,$policy,$j," >> $res
                    srun -n1 --cpus-per-task=64 ./"$precision"gemm_"$library".x $k $k $k >> $res

                done
            done
        done
    done
done