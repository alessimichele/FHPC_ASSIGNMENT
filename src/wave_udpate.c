#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#include <omp.h>

#include "wave_update.h"
#include "io_init.h"

unsigned char* get_next_square(unsigned char *leftop_cell, int size_current_square, int k, unsigned char* grid){
   
        unsigned char* next_square = (unsigned char*)malloc(4*size_current_square + 4); // (size_current_square + 2)**2 - size_current_square**2)
        
        // top edge
        for (int i=0; i< size_current_square + 2; i++){
            next_square[i] = leftop_cell - k - 1 + i;
        }

        // left edge
        for (int i=0; i< size_current_square + 1; i++){
            next_square[size_current_square + 1 + i] = next_square[size_current_square+1] + k*i
        }

        // bottom edge
        for (int i=0; i<size_current_square +1; i++){
            next_square[size_current_square +1 + size_current_square*k] =
        }

    }
}

void smart_static_update(unsigned char* grid, unsigned char* next, int k){


}





    
    // randomly select one cell of the grid to be the source of the wave
    srand(0);
    int rand_cell = rand() % (k*k);

void wave_update(unsigned char* grid, unsigned char* next, int k, int n, int s ){
    
    // randomly select one cell of the grid to be the source of the wave
    srand(0);
    int rand_cell = rand() % (k*k);
    




    // OpenMP implementation upon one process
    for (int step=0; step<n; step++){
        int nthreads;
        #pragma omp parallel
        {
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

                sum = grid[i*k+prev_col] + 
                grid[i*k+next_col] + 
                grid[prev_row*k+j] + 
                grid[next_row*k+j] + 
                grid[prev_row*k+prev_col] + 
                grid[prev_row*k+next_col] + 
                grid[next_row*k+prev_col] + 
                grid[next_row*k+next_col];

                next[i*k+j] = (sum > 765 || sum < 510) ? 0 : 255; 
            }
        }
        }
        // end of parallel region

        unsigned char* tmp;
        tmp = next;
        next = grid;
        grid=tmp;


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
    return;
};
