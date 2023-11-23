#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#include <omp.h>

#include "wave_update.h"
#include "io_init.h"

unsigned int* recoverSquare(int k, int index, int radius) {

    // FREE IDXS MALLOC

  if (radius > (k / 2 -1 )) 
    {
    perror("Radius too large");
    return NULL;
    }
    int tmp = (4*(2*radius - 1) + 4);

    unsigned int* idxs = (unsigned int*)malloc(tmp*sizeof(unsigned int));

    int count=0;
    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {
  
            if (i > -radius && i < radius && j > -radius && j < radius)
                continue; // Skip the inner square
            int row = (index / k + i + k) % k; // "mod k"
            int col = (index % k + j + k) % k;
            

            unsigned int currentIndex = row * k + col;

            idxs[count] = currentIndex;
            count++;
        }
    }
    return idxs;
}

void wave_update(unsigned char* grid, unsigned char* next, int k, int n, int s ){
    
    // randomly select one cell of the grid to be the source of the wave
    srand(0);
    int rand_cell_idx = rand() % (k*k);

    printf("entered wave update\n");

    for(int step=0; step<n; step++){

        printf("Step %d\n", step);

        // questo for itera sui raggi
        for (int radius=1; radius<= (k / 2 + 1 ); radius++){
            printf("Radius %d\n", radius);
            int tmp = (4*(2*radius - 1) + 4);
            unsigned int* idxs = recoverSquare(k, rand_cell_idx, radius);
            
            
            // questo for itera sulle celle del quadrato
            for (int ii=0; ii<tmp;ii++){
                
                printf("Entered inner for\n");
                
                int sum;
                
                int prev_col = (idxs[ii] -1)%(k*k);
                int next_col = (idxs[ii] +1)%(k*k);
                
                sum = grid[prev_col]+
                grid[(prev_col + k)%(k*k)] + 
                grid[(prev_col -k)%(k*k)] +
                grid[(idxs[ii] - k)%(k*k)] + 
                grid[(idxs[ii] + k)%(k*k)] + 
                grid[(next_col - k)%(k*k)] +
                grid[(next_col + k)%(k*k)] +
                grid[next_col];
               
                next[idxs[ii]] = (sum > 765 || sum < 510) ? 0 : 255;  // salvo  per ogni cella del quadrato il suo next state
                }
                printf("check 00\n");
                // sostituisco gli stati aggiornati delle celle del quadrato nella griglia
                for (int ii=0; ii<tmp; ii++){
                    grid[idxs[ii]] = next[idxs[ii]];
                }
                printf("check 1\n");
                free(idxs);
                printf("check 2\n");
        }

/*
        if((step+1)%s==0){
                printf("now  i'm going to write the file\n");

               
                char *file_path = (char*)malloc(30*sizeof(char) + 1);
                strcpy(file_path, "files/wave/");

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
    */
    }

    return;
};







/*


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
    
    
    */