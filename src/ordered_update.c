#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

void ordered_update(unsigned char* grid, int k, int n_steps, int s){
    /*
    evolve the current state of the game of life grid for n_steps using the ordered 
    update algorithm. The grid has size k x k.
    Parameters
    ----------
    grid: pointer to unsigned char, the grid, which is a 1D array, where 
        grid[i*k+j] is the state of the cell in row i and column j
    k: int, size of the grid
    n_steps: int, number of steps to evolve the grid
    s: int, every how many steps a dump of the system is saved on a file
        (0 meaning only at the end)
    */
    for (int step = 0; step < n_steps; step++)
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
        if ((current_step+1) % s == 0)
        {
            char path[45] = "images/evolve_ordered/";
            char name[20];
            snprintf(name, 20, "snapshot_%05d.pgm", current_step+1);
            strcat(path, name);
            //DA RIVEDERE QUANDO AVREMO SCRITTO LA FUNZIONE write_pgm_parallel
            write_pgm_parallel(grid, 255, k, k, path, 0, 1, k);
            
        }
    }
}
