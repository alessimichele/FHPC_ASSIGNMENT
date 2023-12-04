# Exercise 2
This folder contains all the files for the second exercise. The goal of this exercise is to compare performance of three math libraries available on HPC: MKL, OpenBLAS and BLIS. The comparison is performed focusing on the level 3 BLAS function called gemm.

## Files in this folder
- `gemm_parallel.c` is a standard gemm code, where 3 matrices A,B,C are allocated, A and B are filled and the BLAS routine calculates the matrix-matrix product C=A\*B, using parallelization via OpenMP.
- `Makefile` is the makefile used to compile the code.
- `*.sh` are the scripts used to run the code on the cluster.  `$NODE_sbatch.sh` files are used for the scalability over the matrix size, while  `$NODE_cores.sh` files are used for the scalability over the number of cores.
-  `data/ ` directory contains the data produced by the code: it is organized in subdirectories, one for each node used. Each subdirectory contains the data for the scalability over the matrix size, and the scalability over the number of cores.

## Usage
After deciding which node to use and what kind of scalability to test, open the corresponding script and change the location of your BLIS installation. Then, run the script with `sbatch $NODE_$SCALABILITY.sh`. The data will be saved in the corresponding subdirectory of `data/`. The sbatch files already compile the code, so there is no need to do it manually.