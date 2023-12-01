#!/bin/bash 
#SBATCH --partition=THIN
#SBATCH --job-name=lines
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00 
#SBATCH --exclusive
#SBATCH --output=remaining_lines.out

module load mkl
module load openBLAS
export LD_LIBRARY_PATH=/u/dssc/mcarol00/blis/myblis/lib:$LD_LIBRARY_PATH

for i in 1 2 3 4 5
do
    export BLIS_NUM_THREADS=2
    echo "double blis 2 cpus" >> remaining_lines.out
    srun -n1 --cpus-per-task 2 ./dgemm_blis.x 10000 10000 10000 >> remaining_lines.out
    echo "single blis 2 cpus" >> remaining_lines.out
    srun -n1 --cpus-per-task 2 ./fgemm_blis.x 10000 10000 10000 >> remaining_lines.out
    export BLIS_NUM_THREADS=1
    echo "double blis 1 cpu" >> remaining_lines.out
    srun -n1 --cpus-per-task 1 ./dgemm_blis.x 10000 10000 10000 >> remaining_lines.out
    echo "single blis 1 cpu" >> remaining_lines.out
    srun -n1 --cpus-per-task 1 ./fgemm_blis.x 10000 10000 10000 >> remaining_lines.out
    echo "double mkl 2 cpus" >> remaining_lines.out
    srun -n1 --cpus-per-task 2 ./dgemm_mkl.x 10000 10000 10000 >> remaining_lines.out
    echo "single mkl 2 cpus" >> remaining_lines.out
    srun -n1 --cpus-per-task 2 ./fgemm_mkl.x 10000 10000 10000 >> remaining_lines.out
    echo "double mkl 1 cpu" >> remaining_lines.out
    srun -n1 --cpus-per-task 1 ./dgemm_mkl.x 10000 10000 10000 >> remaining_lines.out
    echo "single mkl 1 cpu" >> remaining_lines.out
    srun -n1 --cpus-per-task 1 ./fgemm_mkl.x 10000 10000 10000 >> remaining_lines.out
    echo "double oblas 2 cpus" >> remaining_lines.out
    srun -n1 --cpus-per-task 2 ./dgemm_oblas.x 10000 10000 10000 >> remaining_lines.out
    echo "single oblas 2 cpus" >> remaining_lines.out
    srun -n1 --cpus-per-task 2 ./fgemm_oblas.x 10000 10000 10000 >> remaining_lines.out
    echo "double oblas 1 cpu" >> remaining_lines.out
    srun -n1 --cpus-per-task 1 ./dgemm_oblas.x 10000 10000 10000 >> remaining_lines.out
    echo "single oblas 1 cpu" >> remaining_lines.out
    srun -n1 --cpus-per-task 1 ./fgemm_oblas.x 10000 10000 10000 >> remaining_lines.out
done