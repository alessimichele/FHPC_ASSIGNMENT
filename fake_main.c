#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <mpi.h>



// openMPI/4.1.5/gnu/12.2.1
// salloc --ntasks=3 -N1 -p EPYC --time=00:30:00 
// mpicc -fopenmp main.c -o main
// srun mpirun main  -np 3




int map_even_grid(int rand_cell_idx, int k);
unsigned int* recoverSquare(int k, int index, int radius);
void wave_update_OpenMP(unsigned char* grid, unsigned char* next, int k, int n, int s );
void wave_update_MPI(unsigned char *grid, unsigned char* next, int k, int n, int s, int rank, int size);





/*
----------------------------------------------------------------------------------------------------------------
---------------------------------------------STATIC OpenMP------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
*/





/*
----------------------------------------------------------------------------------------------------------------
---------------------------------------------STATIC MPI------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
*/





/*
----------------------------------------------------------------------------------------------------------------
---------------------------------------------WAVE OpenMP------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
*/


int map_even_grid(int rand_cell_idx, int k){
    /*
    Function that finds the "missing" cell.

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

unsigned int* recoverSquare(int k, int index, int radius) {
    /*
    Function that returns the indexes of the cells that are in the square of radius radius centered in the cell with the given index.

    k: grid size
    index: index of the cell that is the center of the square
    radius: radius of the square
    */

    if (radius > ((k / 2))) {
        perror("Something went wrong. Radius too large");
        return NULL;
    }


    int tmp = (4*(2*radius - 1) + 4); // number of cells in the square of the current radius, i. e. length of indexes

    unsigned int* indxs = (unsigned int*)malloc(tmp*sizeof(unsigned int));

    int count=0;
    #pragma omp parallel
    {
        #pragma omp for
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
    }
    return indxs;
}


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

        /*
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
        */
    } // end of iteration over step
    

    return;
};


/*
----------------------------------------------------------------------------------------------------------------
---------------------------------------------WAVE MPI ------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
*/

void wave_update_MPI(unsigned char *grid, unsigned char* next, int k, int n, int s, int rank, int size){
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

        // broadcast the random cell index
        MPI_Bcast(&rand_cell_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        int thresh; 
        k%2==0 ? (thresh = (k-1)/2) : (thresh =k/2);

        // iterate over the radius
        for (int radius=1; radius<= thresh; radius++){
            
            // rank 0 compute the indexes of the cells in the square          
            int tmp1 = (4*(2*radius - 1) + 4);
            unsigned int* idxs = recoverSquare(k, rand_cell_idx, radius);

            // le prossime 4 righe + il for puo farle il rank 0 e poi broadcastare
            int rem = tmp1 % size;
            int sum = 0;
            int *sendcounts = (int*)malloc(size*sizeof(int));
            int *displs = (int*)malloc(size*sizeof(int));
            

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

            /*
            unsigned int *recvbuf = (unsigned int*)malloc((sendcounts[rank] + 2)*sizeof(unsigned int));
            MPI_Scatterv(idxs, sendcounts, displs, MPI_INT, recvbuf, sendcounts[rank] + 2, MPI_INT, 0, MPI_COMM_WORLD);
            */

            #pragma omp parallel
            {
                #pragma omp for
                // iterate over the given set of cells, according to the process rank
                for (int ii=displs[rank]; ii<sendcounts[rank];ii++){ // non sarà fino a tmp, ma sarà fino a cell_per_proc o qlcs di simile
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

                    /*#pragma omp for
                    for (int ii=displs[rank]; ii<sendcounts[rank]; ii++){ // idem
                        grid[idxs[ii]] = next[idxs[ii]];
                    }*/

                    int *senddispls;
                    if (rank==0){
                        senddispls = (int*)malloc(size*sizeof(int));
                        for (int i=0; i<size; i++){
                            senddispls[i] += sendcounts[i];
                        }
                    }
                    MPI_Bcast(senddispls, size, MPI_INT, 0, MPI_COMM_WORLD);

                    MPI_Barrier(MPI_COMM_WORLD);

                    MPI_Allgatherv(next, sendcounts[rank], MPI_UNSIGNED_CHAR, grid, sendcounts, senddispls, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);


                } // end of parallel region    
                free(idxs);

            free(sendcounts);
            free(displs);

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

        printf("Rank %d has completed step %d.\n", rank, step);

        /*
        if (rank==0){
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
        */

    } // end of iteration over step
   
}


/*
----------------------------------------------------------------------------------------------------------------
--------------------------------------------- FAKE MAIN------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
*/


int main(){
    // initialize the grid
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int k = 100;
    unsigned char *grid = (unsigned char*)malloc(k*k*sizeof(unsigned char));
    unsigned char *next = (unsigned char*)malloc(k*k*sizeof(unsigned char));
    int n = 5;
    int s = 1;


    if (rank==0){
        for (int i=0; i<k*k; i++){
            grid[i] = rand()%2 * 255;
        }
    }
    // broadcast the grid
    MPI_Bcast(grid, k*k, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    



    clock_t start, end;
    double cpu_time_used;

    start = omp_get_wtime();

    wave_update_MPI(grid, next,k, n,  s, rank, size);
    // wave_update_OpenMP(grid, next, k, n, s);
    // static_update_OpenMP(grid, next, k, n, s);

    end = omp_get_wtime();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Time took %f seconds to execute \n", cpu_time_used);
    MPI_Finalize();
    return 0;
}