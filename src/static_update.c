#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#include <omp.h>

#include "static_update.h"
#include "io_init.h"

// current: pointer to the current grid
// next: pointer to the next grid
// k: gride size
// n: number of iterations to be calculated
// s: every s-th iteration a dump of the grid is saved on a file
// rank: rank of the process
// size: number of processes
// rows_per_process: number of rows per process



void static_update_OpenMP(unsigned char *current, unsigned char* next, int k,  int n,  int s){
    
    // OpenMP implementation upon one process
    for (int n=0; n<n; n++){
        int nthreads;
        #pragma omp parallel{
        int id = omp_get_thread_num();

        #pragma omp master
            nthreads = omp_get_num_threads();
        

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


        if((step+1)%s==0){
                printf("now  i'm going to write the file\n");

               
                char *file_path = (char*)malloc(32*sizeof(char) + 1);
                strcpy(file_path, "files/static/");

                char *fname = (char*)malloc(20*sizeof(char) + 1);
                snprintf(fname, 20, "snapshot_%05d.pgm", step+1);
                printf("fname: %s\n", fname);

            
                strcat(file_path, fname);
                // print the file path
                printf("file path: %s\n", file_path);
                printf("address of file_path: %p\n", file_path);

                write_pgm_image((void *)grid, 255, k, k, file_path);

                free(fname);
                free(file_path);
                
            
            }
    }
    }
};





