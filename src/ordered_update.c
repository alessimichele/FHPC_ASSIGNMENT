#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ordered_update.h"
#include "io_init.h"

void ordered_update(unsigned char* grid, int k, int n, int s){
    /*
    evolve the current state of the game of life grid for n using the ordered 
    update algorithm. The grid has size k x k.
    Parameters
    ----------
    grid: pointer to unsigned char, the grid, which is a 1D array, where 
        grid[i*k+j] is the state of the cell in row i and column j
    k: int, size of the grid
    n: int, number of steps to evolve the grid
    s: int, every how many steps a dump of the system is saved on a file
        (0 meaning only at the end)
    */


    for (int step = 0; step < n; step++)
    {   
        for (int i = 0; i < k; i++) // Loop over all rows
        { 
            for (int j = 0; j < k; j++) // Loop over all columns
            {   
                int next_row = (i+1+k)%k;
                int previous_row = (i-1+k)%k;
                int next_column = (j+1+k)%k;
                int previous_column = (j-1+k)%k;
                int n_neigh_255 = grid[previous_row + previous_column] + 
                                  grid[previous_row + next_column] +
                                  grid[next_row + previous_column] + 
                                  grid[next_row + next_column] + 
                                  grid[previous_row + j] +
                                  grid[next_row + j] + 
                                  grid[i*k + previous_column] + 
                                  grid[i*k + next_column];
                // Update cell, if it is alive and has less than 2 or more than 3 neighbours, it dies; if it is dead and has 2 or 3 neighbours, it becomes alive
                grid[i*k+j] = (n_neigh_255 > 765 || n_neigh_255 < 510) ? 0 : 255;
            }
        }
        

            if((step+1)%s==0){
                //printf("now  i'm going to write the file\n");

               
                char *file_path = (char*)malloc(32*sizeof(char) + 1);
                strcpy(file_path, "files/ordered/");

                char *fname = (char*)malloc(20*sizeof(char) + 1);
                snprintf(fname, 20, "snapshot_%05d.pgm", step+1);
                //printf("fname: %s\n", fname);

            
                strcat(file_path, fname);
                // print the file path
                //printf("file path: %s\n", file_path);
                //printf("address of file_path: %p\n", file_path);

                //write_pgm_image((void *)grid, 255, k, k, file_path);

                free(fname);
                free(file_path);
                
            
            }
            
    }
}
