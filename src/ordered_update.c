#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include "ordered_update.h"
#include "io_init.h"

void ordered_update(unsigned char *grid, int k,  int n,  int s, int rank, int size, int my_rows_number){
    if (size == 1){
        ordered_update_OpenMP(grid, k, n, s, rank, size, my_rows_number);
    }else{
        ordered_update_MPI(grid, k, n, s, rank, size, my_rows_number);
    };
}

void ordered_update_OpenMP(unsigned char* grid, int k, int n, int s, int rank, int size, int my_rows_number){
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
   int nthreads;
    #pragma omp parallel
    {
        #pragma omp master
        {
            nthreads = omp_get_num_threads();
        }
    }

    double t_ordered;
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){t_ordered = omp_get_wtime();}

    for (int step = 0; step < n; step++){   
        #pragma omp parallel for ordered
        for (int i = 0; i < k; i++){ // Loop over all rows
            #pragma omp ordered
            for (int j = 0; j < k; j++){ // Loop over all columns
                int next_row = (i+1+k)%k;
                int previous_row = (i-1+k)%k;
                int next_column = (j+1+k)%k;
                int previous_column = (j-1+k)%k;
                int sum = grid[previous_row + previous_column] + 
                    grid[previous_row + next_column] +
                    grid[next_row + previous_column] + 
                    grid[next_row + next_column] + 
                    grid[previous_row + j] +
                    grid[next_row + j] + 
                    grid[i*k + previous_column] + 
                    grid[i*k + next_column];
                // Update cell, if it is alive and has less than 2 or more than 3 neighbours, it dies; if it is dead and has 2 or 3 neighbours, it becomes alive
                grid[i*k+j] = (sum > 765 || sum < 510) ? 0 : 255;
            }
        }
        if((step+1)%s==0){
            char *file_path = (char*)malloc(32*sizeof(char) + 1);
            strcpy(file_path, "files/ordered/");

            char *fname = (char*)malloc(20*sizeof(char) + 1);
            snprintf(fname, 20, "snapshot_%05d.pgm", step+1);
            strcat(file_path, fname);

            parallel_write(grid, 255, file_path, k, my_rows_number, rank, size, MPI_COMM_WORLD);

            free(fname);
            free(file_path);
        }    
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){t_ordered = omp_get_wtime() - t_ordered;
        printf("r,%d,%d,%d,%lf\n", size, nthreads, k, t_ordered);
    }
    return;
}
void ordered_update_MPI(unsigned char* grid, int k, int n, int s, int rank, int size, int my_rows_number){
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
   int nthreads;
    #pragma omp parallel
    {
        #pragma omp master
        {
            nthreads = omp_get_num_threads();
        }
    }

    double t_ordered;
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){t_ordered = omp_get_wtime();}

 
    unsigned char* previous_row = (unsigned char*)malloc(k*sizeof(unsigned char));
    unsigned char* next_row = (unsigned char*)malloc(k*sizeof(unsigned char));
    int next_tag = 1000000;
    for (int step=0; step<n;step++){
        MPI_Status status;
        MPI_Request request;
        //send operations have tag = step*size + rank
        //receive operations have tag = step*size + rank +- 1
        if (rank==size-1){
            MPI_Send(grid + k*(my_rows_number-1), k, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
               
        }
        if (rank!=0){
            //send my first row to rank-1
            MPI_Isend(grid, k, MPI_UNSIGNED_CHAR, rank-1, next_tag+ step*size +(rank), MPI_COMM_WORLD, &request);
            //receive my last row from rank+1
            MPI_Recv(next_row, k, MPI_UNSIGNED_CHAR, (rank + 1)%size,next_tag+ step*size+((rank+1)%size), MPI_COMM_WORLD, &status);
            //receive previous row from rank-1
            MPI_Recv(previous_row, k, MPI_UNSIGNED_CHAR, rank - 1, step*size+(rank-1), MPI_COMM_WORLD, &status);
        
        }else{
            MPI_Recv(previous_row, k, MPI_UNSIGNED_CHAR, size-1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(next_row, k, MPI_UNSIGNED_CHAR, 1,next_tag+ step*size+(rank+1), MPI_COMM_WORLD, &status);
        }

        //update the first row
        for (int j=0; j<k; j++){
            int sum;
            int prev_col = (j - 1 + k)%k;
            int next_col = (j + 1 + k)%k;
            int next_row = 1;
            sum = previous_row[prev_col] + 
            previous_row[j] +
            previous_row[next_col] + 
            grid[next_row*k+j] + 
            grid[next_row*k+prev_col] + 
            grid[next_row*k+next_col] +
            grid[0*k+prev_col] +
            grid[0*k+next_col];
            grid[j] = (sum > 765 || sum < 510) ? 0 : 255; 
        }
    

        for (int i=1; i<my_rows_number-1; i++){
            for (int j=0; j<k; j++){
                int sum;
                int prev_col = (j - 1 + k)%k;
                int next_col = (j + 1 + k)%k;
                int prev_row = i-1;
                int next_row = i+1;
                sum = grid[prev_row*k+prev_col] + 
                grid[prev_row*k+j] +
                grid[prev_row*k+next_col] +
                grid[next_row*k+prev_col] + 
                grid[next_row*k+next_col] +
                grid[next_row*k+j] +
                grid[i*k+prev_col] +
                grid[i*k+next_col];
                grid[i*k+j] = (sum > 765 || sum < 510) ? 0 : 255; 
            }
        }
        //update last row
        for (int j=0; j<k; j++){
            int sum;
            int prev_col = (j - 1 + k)%k;
            int next_col = (j + 1 + k)%k;
            int prev_row = my_rows_number-2;
            sum = grid[prev_row*k+prev_col] + 
            grid[prev_row*k+j] +
            grid[prev_row*k+next_col] +
            grid[(my_rows_number-1)*k+prev_col] + 
            grid[(my_rows_number-1)*k+next_col] +
            next_row[prev_col] +
            next_row[next_col] +
            next_row[j];
            grid[(my_rows_number-1)*k+j] = (sum > 765 || sum < 510) ? 0 : 255; 
        }
    
        if (rank!=size-1){
            //send my last row to rank+1
            MPI_Isend(grid + k*(my_rows_number-1), k, MPI_UNSIGNED_CHAR, rank+1, step*size +(rank), MPI_COMM_WORLD, &request);
        }

        if (rank == 0){
            MPI_Send(grid, k, MPI_UNSIGNED_CHAR, size-1,next_tag+ step*size +(rank), MPI_COMM_WORLD);
        }
        
        if ((step+1)%s==0){
            char *file_path = (char*)malloc(32*sizeof(char) + 1);
            strcpy(file_path, "files/ordered/");

            char *fname = (char*)malloc(20*sizeof(char) + 1);
            snprintf(fname, 20, "snapshot_%05d.pgm", step+1);
            strcat(file_path, fname);

            parallel_write(grid, 255, file_path, k, my_rows_number, rank, size, MPI_COMM_WORLD);
    
            free(fname);
            free(file_path);  
        }
       
    }
    free(previous_row);
    free(next_row);

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){t_ordered = omp_get_wtime() - t_ordered;
        printf("r,%d,%d,%d,%lf\n", size, nthreads, k, t_ordered);
    }
    return;
}