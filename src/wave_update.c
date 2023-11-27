#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>

#include "wave_update.h"
#include "io_init.h"

/*
----------------------------------------------------------------------------------------------------------------
---------------------------------------------WRAPPER------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
*/

void wave_update(unsigned char *grid, unsigned char* next, int k,  int n,  int s, int size, int rank){
    /*
    Wrapper function for the wave update

    grid: pointer to the grid grid
    next: pointer to the next grid
    k: gride size
    n: number of iterations to be calculated
    s: every s-th iteration a dump of the grid is saved on a file
    size: number of processes
    rank: rank of the process
    */

    if (size == 1){
        wave_update_OpenMP(grid, next, k, n, s);    
    }else{
        wave_update_MPI(grid, next, k, n, s, size, rank);
    }
    return;
}


/*
----------------------------------------------------------------------------------------------------------------
---------------------------------------------AUXILIARY FUNCTIONS------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
*/



int map_even_grid(int rand_cell_idx, int k){
    /*
    Function that finds the cell in the intersection of the missing row and the missing col in case k even.

    rand_cell_idx: index of the cell that is the source of the wave
    k: grid size
    */
    if (k%2 !=0 ){
        perror("Something went wrong, this function should be called only when k is even.\n");
    }
    if (rand_cell_idx%k < 2){
        return (rand_cell_idx+ ((k/2) *k))%(k*k) + k/2;
    }else{
        return (rand_cell_idx+ ((k/2) *k))%(k*k) - k/2;
    }
}

unsigned int* recoverSquare(int k, int index, int radius){
    if (radius > ((k / 2))) {
        perror("Something went wrong. Radius too large");
        return NULL;
    }

    int tmp = (4*(2*radius - 1) + 4);

    unsigned int* indxs = (unsigned int*)malloc(tmp*sizeof(unsigned int));

    int a = index/k; 
    int b = index%k;

    // top-left corner
    int t_l = k * ((a - radius + k) % k) + (b - radius + k) % k;

    // top-right corner
    int t_r = k * ((a - radius + k) % k) + (b + radius) % k;

    // bottom-left corner
    int b_l = k * ((a + radius) % k) + (b - radius + k) % k;

    // bottom-right corner
    int b_r = k * ((a + radius) % k) + (b + radius) % k;

    int l_edge = tmp/4;

    // top edge
    for (int ii=0; ii<l_edge; ii++) {
        int row = t_l/k;
        indxs[ii] = k*row + (t_l + ii) % k;
    }

    // right edge
    for (int ii=0; ii<l_edge; ii++) {
        int col = t_r%k;
        indxs[ii + l_edge] = (t_r + ii*k )%(k*k);
    }

    // bottom edge
    for (int ii=0; ii<l_edge; ii++) {
        int row = b_r/k;
        indxs[ii + 2*l_edge] = k*row + (b_r - ii + k) % k;
    }

    // left edge
    for (int ii=0; ii<l_edge; ii++) {
        int col = b_l%k;
        indxs[ii + 3*l_edge] = (b_l - ii*k +(k*k))%(k*k);
    }

    // sort the idxs array
    for (int i=0; i<tmp; i++){
        for (int j=i+1; j<tmp; j++){
            if (indxs[i] > indxs[j]){
                int tmp = indxs[i];
                indxs[i] = indxs[j];
                indxs[j] = tmp;
            }
        }
    }

    return indxs;
}


/*
----------------------------------------------------------------------------------------------------------------
---------------------------------------------OpenMP UPDATE------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
*/

void wave_update_OpenMP(unsigned char* grid, unsigned char* next, int k, int n, int s ){
    /*
    Function that updates the grid according to the wave rule.

    grid: pointer to the grid grid
    next: pointer to the next grid
    k: gride size
    n: number of iterations to be calculated
    s: every s-th iteration a dump of the grid is saved on a file
    */

    //printf("entered wave update\n");

    // iterate over the step
    for(int step=0; step<n; step++){

        // randomly select one cell of the grid to be the source of the wave
        // generate a seed depending on millisec
        struct timeval time;
        gettimeofday(&time, NULL);
        int millis = ((time.tv_sec * 1000) + (time.tv_usec / 1000))%1000;
        srand(millis);
        int rand_cell_idx = rand() % (k*k); // rand_cell_idx is the index of the cell that will be the source of the wave

        int thresh; // threshold for the maximum possible radius, depends on if k is even or odd
        k%2==0 ? (thresh = (k-1)/2) : (thresh =k/2);

        // iterate over the radius
        for (int radius=1; radius<= thresh; radius++){
            
            int tmp1 = (4*(2*radius - 1) + 4); // number of cells in the square of the current radius
            unsigned int* idxs = recoverSquare(k, rand_cell_idx, radius);

            
            if (idxs==NULL){
                // go to next iteration of the outer outer loop (the one with step)
                // if code works properly, this should never happen
                printf("idxs is NULL... something is going wrong, check it please.\n");
                continue;
            }

            #pragma omp parallel
            {
                #pragma omp for
                // iterate over the cells in the square
                for (int ii=0; ii<tmp1;ii++){
                    int prev_col = (idxs[ii] -1 +(k*k))%(k*k);
                    int next_col = (idxs[ii] +1 +(k*k))%(k*k);

                    int sum=0;
                    sum += grid[prev_col]+
                    grid[(prev_col + k +(k*k))%(k*k)] +
                    grid[(prev_col -k +(k*k))%(k*k)] +
                    grid[(idxs[ii] - k +(k*k))%(k*k)] + 
                    grid[(idxs[ii] + k +(k*k))%(k*k)] + 
                    grid[(next_col - k +(k*k))%(k*k)] +
                    grid[(next_col + k +(k*k))%(k*k)] +
                    grid[next_col];
                
                    next[idxs[ii]] = (sum > 765 || sum < 510) ? 0 : 255;  // salvo  per ogni cella del quadrato il suo next state
                
                    } // end of iteration over cells in the given square

                    #pragma omp for
                    // sostituisco gli stati aggiornati delle celle del quadrato nella griglia
                    for (int ii=0; ii<tmp1; ii++){
                        grid[idxs[ii]] = next[idxs[ii]];
                    }
                } // end of parallel region    
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

            #pragma omp parallel
            {
            // update last row and col
            #pragma omp for
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
            #pragma omp for
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
        
            #pragma omp for
            for (int i=0; i<k; i++){
                grid[i*k + col] = next[i*k + col];
            }
            #pragma omp for
            for (int j=0; j<k; j++){
                if (j!= col){
                    grid[row*k + j] = next[row*k + j]; 
                }
            }

            } // end of the parallel region 
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
    } // end of iteration over step

    return;
};


/*
----------------------------------------------------------------------------------------------------------------
---------------------------------------------MPI UPDATE------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
*/

void wave_update_MPI(unsigned char *grid, unsigned char* next, int k, int n, int s){
    /*
    Function that updates the grid according to the wave rule, using MPI.

    grid: pointer to the grid grid
    next: pointer to the next grid
    k: gride size
    n: number of iterations to be calculated
    s: every s-th iteration a dump of the grid is saved on a file
    */

    // iterate over the step
    for (int step=0; step<n; step++){

        if (rank==0 && step%10==0){
            printf("----------------Step %d--------------\n", step);       
        }
        
        // rank 0 compute the random cell index
        // randomly select one cell of the grid to be the source of the wave
        int rand_cell_idx;
        if (rank == 0){
            struct timeval time;
            gettimeofday(&time, NULL);
            int millis = ((time.tv_sec * 1000) + (time.tv_usec / 1000))%1000;
            srand(millis);
            rand_cell_idx = rand() % (k*k); 
        }

        // broadcast the random cell index to all the processes
        MPI_Bcast(&rand_cell_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        int thresh; 
        k%2==0 ? (thresh = (k-1)/2) : (thresh =k/2);

        MPI_Barrier(MPI_COMM_WORLD); // wait for all processes to get the random cell index

        // iterate over the radius
        for (int radius=1; radius<= thresh; radius++){
            
                      
            int tmp1 = (4*(2*radius - 1) + 4);
            unsigned int* idxs = recoverSquare(k, rand_cell_idx, radius);

            int rem = tmp1 % size;
            int sum = 0;
            int *sendcounts = (int*)calloc(size*sizeof(int), sizeof(int));
            int *displs = (int*)calloc(size*sizeof(int), sizeof(int));
            

            // calculate sendcounts and displs
            for (int i=0; i<size; i++){
                sendcounts[i] = tmp1/size;
                if (rem > 0){
                    sendcounts[i]++;
                    rem--;
                }
                displs[i] = sum;
                sum += sendcounts[i];
            }


         
            #pragma omp parallel for
            // iterate over the given set of cells, according to the process rank
            for (int ii=displs[rank]; ii<sendcounts[rank];ii++){ 
                int prev_col = (idxs[ii] -1 +(k*k))%(k*k);
                int next_col = (idxs[ii] +1 +(k*k))%(k*k);

                int sum=0;
                sum += grid[prev_col]+
                grid[(prev_col + k +(k*k))%(k*k)] +
                grid[(prev_col -k +(k*k))%(k*k)] +
                grid[(idxs[ii] - k +(k*k))%(k*k)] + 
                grid[(idxs[ii] + k +(k*k))%(k*k)] + 
                grid[(next_col - k +(k*k))%(k*k)] +
                grid[(next_col + k +(k*k))%(k*k)] +
                grid[next_col];
            
                next[idxs[ii]] = (sum > 765 || sum < 510) ? 0 : 255;  
            
                } // end of iteration over cells

            MPI_Barrier(MPI_COMM_WORLD);

            int offset;
            if (rank!=size-1){
                offset = idxs[displs[rank+1]] - idxs[displs[rank]] - 1;
            }else{
                offset =  idxs[tmp1-1] - idxs[displs[rank]] + 1;
            }

            int *recvcounts = (int*)calloc(size*sizeof(int), sizeof(int));
            for (int i=0; i<size-1; i++){
                recvcounts[i] = idxs[displs[i+1]] - idxs[displs[i]];
            }
            recvcounts[size-1] = idxs[tmp1-1] - idxs[displs[size-1]] + 1;

            int *recvdispls = (int*)calloc(size*sizeof(int), sizeof(int));
            for (int i=0; i<size; i++){
                recvdispls[i] = idxs[displs[i]];
            }

            MPI_Allgatherv(next + idxs[displs[rank]], offset, MPI_UNSIGNED_CHAR, grid, recvcounts, recvdispls, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD); 

            free(idxs);
            free(sendcounts);
            free(displs);
            free(recvcounts);
            free(recvdispls);
        } // end of iteration over radius

        // FLAG
        if (k%2==0){
            if (rank==0 || rank==size-1){
                int crucial_point = map_even_grid(rand_cell_idx, k);
                int row = crucial_point / k;
                int col = crucial_point % k;

                if (rank==size-1){ // last process updates the missing row
                    #pragma omp parallel for
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

                    #pragma omp parallel for
                    for (int i=0; i<k; i++){
                        grid[i*k + col] = next[i*k + col];
                    }

                }else if (rank==0){ // first process updates the missing col
                    #pragma omp parallel for
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
                    #pragma omp parallel for
                    for (int j=0; j<k; j++){
                        if (j!= col){
                            grid[row*k + j] = next[row*k + j]; 
                        }
                    } 
                }

                if(size==1){
                    perror("Something went wrong. 'wave_update_MPI' should be called only when size is greater than 1.\n");
                    return;
                }

                if(rank==size-1){ // last process sends the updated row to the first process
                    MPI_Send(grid + row, k, MPI_UNSIGNED_CHAR, 0, step, MPI_COMM_WORLD);
                }else if(rank==0){ // first process receives the updated row from the last process
                    MPI_Recv(grid + row, k, MPI_UNSIGNED_CHAR, size-1, step, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

            MPI_Bcast(grid, k*k, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD); // wait for all processes to get the updated grid
        } // end of the flag

        printf("step %d completed\n", step);

        /*if (rank==0){
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
        }*/

    } // end of iteration over step

     printf("Rank %d has completed all the steps.\n", rank);
    return;
}
