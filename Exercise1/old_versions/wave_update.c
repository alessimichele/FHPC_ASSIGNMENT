#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#include <omp.h>

#include "wave_update.h"
#include "io_init.h"


int map_even_grid(int rand_cell_idx, int k){
    if (k%2 !=0 ){
        perror("Something went wrong, this function should be called only when k is even.\n");
    }
    if (rand_cell_idx%k < 2){
        return (rand_cell_idx+ ((k/2) *k))%(k*k) + k/2;
    }else{
        return (rand_cell_idx+ ((k/2) *k))%(k*k) - k/2;
    }
}

unsigned int* recoverSquare(int k, int index, int radius) {
  if (radius > ((k / 2))) 
    {
    perror("Something went wrong. Radius too large");
    return NULL;
    }

    int tmp = (4*(2*radius - 1) + 4);

    unsigned int* indxs = (unsigned int*)malloc(tmp*sizeof(unsigned int));

    int count=0;
    for (int i = -radius; i <= radius; ++i) {
        for (int j = -radius; j <= radius; ++j) {
  
            if (i > -radius && i < radius && j > -radius && j < radius)
                continue; // Skip the inner square
            int row = (index / k + i + k) % k; // "mod k"
            int col = (index % k + j + k) % k;
            
            unsigned int currentIndex = row * k + col;

            indxs[count] = currentIndex;
            count++;
        }
    }
    return indxs;
}

void wave_update(unsigned char* grid, unsigned char* next, int k, int n, int s ){
    
    // randomly select one cell of the grid to be the source of the wave
    srand(0);
    int rand_cell_idx = rand() % (k*k);

    //printf("entered wave update\n");

    // iterate over the step
    for(int step=0; step<n; step++){

        int thresh;
        k%2==0 ? (thresh = (k-1)/2) : (thresh =k/2);

        // questo for itera sui raggi
        for (int radius=1; radius<= thresh; radius++){
            // printf("Radius %d\n", radius);
            int tmp1 = (4*(2*radius - 1) + 4);
            unsigned int* idxs = recoverSquare(k, rand_cell_idx, radius);

            
            if (idxs==NULL){
                // go to next iteration of the outer outer loop (the one with step)
                // if code works properly, this should never happen
                printf("idxs is NULL... something is going wrong, check it please.\n");
                continue;
            }
            // questo for itera sulle celle del quadrato
            for (int ii=0; ii<tmp1;ii++){
                int prev_col = (idxs[ii] -1+(k*k))%(k*k);
                int next_col = (idxs[ii] +1+(k*k))%(k*k);

                int sum=0;
                sum += grid[prev_col]+
                grid[(prev_col + k+(k*k))%(k*k)] +
                grid[(prev_col -k+(k*k))%(k*k)] +
                grid[(idxs[ii] - k+(k*k))%(k*k)] + 
                grid[(idxs[ii] + k+(k*k))%(k*k)] + 
                grid[(next_col - k+(k*k))%(k*k)] +
                grid[(next_col + k+(k*k))%(k*k)] +
                grid[next_col];
               
                next[idxs[ii]] = (sum > 765 || sum < 510) ? 0 : 255;  // salvo  per ogni cella del quadrato il suo next state
             
                } // end of iteration over cells in the given square
               
                // sostituisco gli stati aggiornati delle celle del quadrato nella griglia
                for (int ii=0; ii<tmp1; ii++){
                    grid[idxs[ii]] = next[idxs[ii]];
                }
                
                free(idxs); 
        } // end of iteration over radius
        

        // FLAG
        if (k%2==0){
            // enter here only if k is even
            // and finish to update the grid

            int crucial_point = map_even_grid(rand_cell_idx, k); // intersection of the missing col and row
            
            // row index and col index of the cell that is the intersection of the non-updated row and col at the end
            int row = crucial_point / k;
            int col = crucial_point % k;
           
            // update last row and col
            for (int i = 0; i<k;i++){
                int j = col;
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
            for (int j = 0; j<k;j++){
                if (j != col){ 
                    int i = row;
                    
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
        

            for (int i=0; i<k; i++){
                grid[i*k + col] = next[i*k + col];
            }
            for (int j=0; j<k; j++){
                if (j!= col){
                    grid[row*k + j] = next[row*k + j]; 
                }
            }   
        } // end of the flag

        printf("step %d completed\n", step);


        if((step+1)%s==0){
                //printf("now  i'm going to write the file\n");
               
                char *file_path = (char*)malloc(29*sizeof(char) + 1);
                strcpy(file_path, "files/wave/");

                char *fname = (char*)malloc(20*sizeof(char) + 1);
                snprintf(fname, 20, "snapshot_%05d.pgm", step+1);
                //printf("fname: %s\n", fname);

                strcat(file_path, fname);
                // print the file path
                //printf("file path: %s\n", file_path);
                //printf("address of file_path: %p\n", file_path);

                write_pgm_image((void *)grid, 255, k, k, file_path);

                free(fname);
                free(file_path);
            }
            
    
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