#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "io_init.h"
#include "ordered_update.h"
#include "static_update.h"
#include "wave_update.h"

#define INIT 1
#define RUN  2
#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9);

#define K_DFLT 1000

#define ORDERED 0
#define STATIC  1
#define WAVE 2
#define ORDERED_FINITE 3

int   action = INIT;
int   k      = K_DFLT;
int   e      = ORDERED;
int   n      = 100;
int   s      = 0;
char *fname  = NULL;
int maxval = 65535; // maybe

int main ( int argc, char **argv ){

    int action = 0;
    char *optstring = "irk:e:f:n:s:";

    int c;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch(c) {
        
        case 'i':
            action = INIT; break;
            
        case 'r':
            action = RUN; break;
            
        case 'k':
            k = atoi(optarg); break;

        case 'e':
            e = atoi(optarg); break;

        case 'f':
            fname = (char*)malloc( sizeof(optarg)+1); //fname now is an "array of char" variable (a string)
            sprintf(fname, "%s", optarg );
            break;

        case 'n':
            n = atoi(optarg); break;

        case 's':
            s = atoi(optarg); break;

        default :
            printf("argument -%c not known\n", c ); break;
        }
    }

    if (s==0){
        s = n; // print last step
    }

  
    

    if(fname == NULL){
        fprintf(stderr, "Filename is not provided. Please provide a filename with -f option. This will be the name of the file containing the initial grid.\n");
        return 1;
    }

    char *file_path = (char*)malloc(sizeof(char)*(strlen(fname)+strlen("./files/init/") + 1));

    strcpy(file_path, "./files/init/");
    strcat(file_path, fname);



    int rank, size;
    MPI_Init( NULL, NULL );
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank( comm,&rank );
    MPI_Comm_size( comm,&size );

    int my_rows_number = (rank<(k%size)) ? k/size + 1 : k/size;

    int nthreads;
    #pragma omp parallel
    {
        #pragma omp master
        {
            nthreads = omp_get_num_threads();
        }
    }


    if (action == INIT){

        if(rank==0)printf("mode,size,nthreads,k,time\n");
        double t_init;
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==0)t_init = omp_get_wtime();
        init_parallel(file_path, k, rank, size, my_rows_number);
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==0){t_init = omp_get_wtime() - t_init;
            printf("i,%d,%d,%d,%lf\n", size, nthreads, k, t_init);
        }

    }else{  

        if (action != RUN) {
            if(rank==0)fprintf(stderr, "Please provide an action. Either -i for initialization or -r for running the simulation.\n");
            return 1; 
        }
    
        unsigned char *partial_grid;
        int maxval, xsize, ysize;
        if(rank==0)printf("mode,size,nthreads,k,time\n");
        double t_read;
        if (rank==0)t_read = omp_get_wtime();
        parallel_read(&partial_grid, &maxval, &xsize, &ysize, file_path, my_rows_number, rank, size, comm); // !!!!!!
        if(rank == 0){t_read = omp_get_wtime() - t_read;
            printf("re,%d,%d,%d,%lf\n", size, nthreads, k, t_read);
        }


        if (e == ORDERED){

            //if(rank==0)printf("Run in order mode.\n");
            ordered_update(partial_grid, k, n, s, rank, size, my_rows_number);
            //if(rank==0)printf("Done!\n");

        }else if (e == STATIC){

            //if(rank==0)printf("Run in static mode.\n");
            unsigned char* next = (unsigned char*)malloc(k*my_rows_number*sizeof(unsigned char));
            static_update(partial_grid, next, k, n, s, rank, size, my_rows_number);
            free(next);
            //if(rank==0)printf("Done!\n");

        }else if (e == WAVE){

            //if(rank==0){printf("Run in wave mode.\n");}
            unsigned char* grid = (unsigned char*)malloc(k*k*sizeof(unsigned char));
            int *recvcountsgather = (int*)malloc(size*sizeof(int));
            int *displsgather = (int*)malloc(size*sizeof(int));
            for(int i=0; i<size; i++){
                recvcountsgather[i] = (i<(k%size)) ? k/size + 1 : k/size;
                recvcountsgather[i] = recvcountsgather[i]*k;
                displsgather[i] = (i==0) ? 0 : displsgather[i-1] + recvcountsgather[i-1];
            }
            MPI_Barrier(MPI_COMM_WORLD);
            
            MPI_Allgatherv(partial_grid, k*my_rows_number, MPI_UNSIGNED_CHAR, grid, recvcountsgather, displsgather, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
            
            unsigned char* next = (unsigned char*)malloc(k*k*sizeof(unsigned char));
            wave_update(grid, next, k, n, s, rank, size, my_rows_number);
            
            free(next);
            free(grid);
            free(recvcountsgather);
            free(displsgather);
            //if(rank==0){printf("Done!\n");}

        }else{
            fprintf(stderr, "Please provide a valid execution mode. Either -e 0 for ordered, -e 1 for static or -e 2 for wave.\n");
        }
        
        free(partial_grid);
    
    }

    free(file_path);
    free(fname);
    MPI_Finalize();
    return 0;
}

