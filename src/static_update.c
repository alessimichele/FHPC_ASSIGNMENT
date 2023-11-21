#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "static_update.h"

// current: pointer to the current grid
// next: pointer to the next grid
// k: gride size
// n: number of iterations to be calculated
// s: every s-th iteration a dump of the grid is saved on a file
// rank: rank of the process
// size: number of processes
// rows_per_process: number of rows per process

void static_update(unsigned char *current, unsigned char* next, int k,  int n,  int s, int rank, int size, int rows_per_process){

    if (size==1){
        static_update_OpenMP(current, next, k, n, s);
    } else {
        static_update_MPI(current, next, k, n, s, rank, size, rows_per_process);
    }
};

void static_update_OpenMP(unsigned char *current, unsigned char* next, int k,  int n,  int s){
    
    // OpenMP implementation upon one process
    for (int n=0; n<n; n++){
        int nthreads;
        #pragma omp parallel
        int id = omp_get_thread_num();

        #pragma omp master{
            nthreads = omp_get_num_threads();
        }

        #pragma omp for
        for (int i=0; i<k; i++){
            for (int j=0; j<k; j++){
                int sum;
                int prev_col = (j - 1 + k)%k;
                int next_col = (j + 1 + k)%k;
                int prev_row = (i - 1 + k)%k;
                int next_row = (i + 1 + k)%k;

                sum = current[i*k+prev_col] + 
                current[i*k+next_col] + 
                current[prev_row*k+j] + 
                current[next_row*k+j] + 
                current[prev_row*k+prev_col] + 
                current[prev_row*k+next_col] + 
                current[next_row*k+prev_col] + 
                current[next_row*k+next_col];

                next[i*k+j] = (sum > 765 || sum < 510) ? 0 : 255; 
            }
        }

        unsigned char* tmp;
        tmp = next;
        next = current;
        current=tmp;


        if ((current_step+1) % s == 0)
        {
            char path[45] = "files/static_update/";
            char name[20];
            snprintf(name, 20, "snapshot_%05d.pgm", current_step+1);
            strcat(path, name);
            // DA RIVEDERE QUANDO AVREMO SCRITTO LA FUNZIONE write_pgm_parallel
            write_pgm_parallel(grid, 255, k, k, path, 0, 1, k);
            
        }
    }

};

void static_update_MPI(unsigned char *current, unsigned char* next, int k,  int n,  int s, int rank, int size, int rows_per_process){


};




