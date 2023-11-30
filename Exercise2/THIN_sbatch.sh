#!/bin/bash 
#SBATCH --partition=THIN
#SBATCH --job-name=blis
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=blis.out

module load architecture/Intel
module load mkl
module load openBLAS/0.3.23-omp
export LD_LIBRARY_PATH=/u/dssc/mcarol00/blis/myblis/lib:$LD_LIBRARY_PATH

export OMP_PLACES=cores
export OMP_NUM_THREADS=12
export BLIS_NUM_THREADS=12

name=complete_scalability_study_THIN_blis
res=THIN/data/$name.csv
echo "library,precision,k,policy,iteration,time,GFLOPS" > $res

srun -n1 make cpu # Now I have all the needed executables.

for library in 'blis' #'mkl' 'blis' 'oblas'
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
                    srun -n1 --cpus-per-task=12 ./"$precision"gemm_"$library".x $k $k $k >> $res

                done
            done
        done
    done
done