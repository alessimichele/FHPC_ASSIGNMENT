#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "static_update.h"
#include "io_init.h"

// grid: pointer to the grid grid
// next: pointer to the next grid
// k: gride size
// n: number of iterations to be calculated
// s: every s-th iteration a dump of the grid is saved on a file
// rank: rank of the process
// size: number of processes
// my_rows_number: number of rows per process

//let's make a wrapper for the static update function
void static_update(unsigned char *grid, unsigned char* next, int k,  int n,  int s, int rank, int size, int my_rows_number){
    if (size == 1){
        static_update_OpenMP(grid, next, k, n, s, rank, size, my_rows_number);
    }else{
        static_update_MPI(grid, next, k, n, s, rank, size, my_rows_number);
    };
}

void static_update_OpenMP(unsigned char *grid, unsigned char* next, int k,  int n,  int s, int rank, int size, int my_rows_number){
    
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

                parallel_write(grid, 255, file_path, k, k, 0, 1, MPI_COMM_WORLD);

                free(fname);
                free(file_path);
            }
    }
};


void static_update_MPI(unsigned char* grid, unsigned char* next, int k, int n, int s, int rank, int size, int my_rows_number){
    
    MPI_Request request[4];  
    unsigned char* previous_row = (unsigned char*)malloc(k*sizeof(unsigned char));
    unsigned char* next_row = (unsigned char*)malloc(k*sizeof(unsigned char));

    for (int step=0; step<n; step++){    
        
        
        //non blocking 
        //send the last row, tag = step*size+rank
        MPI_Isend(grid+(my_rows_number-1)*k, k, MPI_UNSIGNED_CHAR, (rank+1)%size, step*size+rank, MPI_COMM_WORLD, &request[0]); 
        //send the first row
        MPI_Isend(grid, k, MPI_UNSIGNED_CHAR, (rank+size-1)%size, step*size+rank, MPI_COMM_WORLD, &request[1]); 

        //receive the last row
        MPI_Irecv(previous_row, k, MPI_UNSIGNED_CHAR, (rank+size-1)%size, step*size+(rank+size-1)%size, MPI_COMM_WORLD, &request[2]);
        //receive the first row
        MPI_Irecv(next_row, k, MPI_UNSIGNED_CHAR, (rank+1)%size, step*size+(rank+1)%size, MPI_COMM_WORLD, &request[3]);

        MPI_Waitall(4, request, MPI_STATUS_IGNORE);
      
        //update first row
        #pragma omp parallel
        {    
            #pragma omp for   
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

                next[j] = (sum > 765 || sum < 510) ? 0 : 255; 
            }
            
            //update last row
            #pragma omp for
               
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
                next[(my_rows_number-1)*k+j] = (sum > 765 || sum < 510) ? 0 : 255; 
            }
            
            //update the rest of the rows
            #pragma omp parallel for
            
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

                    next[i*k+j] = (sum > 765 || sum < 510) ? 0 : 255; 
                }
            }
            
        }
        //let's syncronize the processes
        MPI_Barrier(MPI_COMM_WORLD);
        //let's switch the pointers
        unsigned char* tmp;
        tmp = next;
        next = grid;
        grid = tmp;
        printf("Rank %d has completed step %d.\n", rank, step);    

        if(rank==0)printf("s: %d\n", s);
        
        if ((step+1)%s==0){
            char *file_path = (char*)malloc(32*sizeof(char) + 1);
            strcpy(file_path, "files/static/");

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
    
    return;
};
