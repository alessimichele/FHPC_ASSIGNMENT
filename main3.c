#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#include "io_init.h"
#include "ordered_update.h"
#include "static_update.h"
#include "wave_update.h"

#define INIT 1
#define RUN  2

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
        s = n;
    }

    printf("k: %d\nn: %d\ns: %d\n", k, n, s);

    printf("fname: %s\n", fname);

    if(fname == NULL){
        fprintf(stderr, "Filename is not provided. Please provide a filename with -f option. This will be the name of the file containing the initial grid.\n");
        return 1;
    }

    char *file_path = (char*)malloc(sizeof(char)*(strlen(fname)+strlen("./files/init/") + 1));

    strcpy(file_path, "./files/init/");
    strcat(file_path, fname);

    printf("file_path: %s\n", file_path);

    int rank, size;
    MPI_Init( NULL, NULL );
    MPI_Comm_rank( MPI_COMM_WORLD,&rank );
    MPI_Comm_size( MPI_COMM_WORLD,&size );

    int my_rows_number = (rank<(k%size)) ? k/size+1 : k/size;

    if (action == INIT){

        if(rank==0)printf("Initializing...\n");
        init_parallel(file_path, k, rank, size, my_rows_number);
        if(rank==0)printf("Initialization done!\n");

    }else{  

        if (action != RUN) {
            if(rank==0)fprintf(stderr, "Please provide an action. Either -i for initialization or -r for running the simulation.\n");
            return 1; 
        }

        if (file_path==NULL){ // questo penso si possa togliere... il punto sarebbe checkare che il file esista
            if(rank==0)fprintf(stderr,"Initial grid not found. Please run with -i for initialization.");
            return 1;
        }
    
        unsigned char *partial_grid;
        int maxval, xsize, ysize;
        parallel_read(&partial_grid, &maxval, &xsize, &ysize, file_path, my_rows_number, rank, size, MPI_COMM_WORLD); // !!!!!!
        
        if(rank==0)printf("The initial grid has been read.\n");

        if (e == ORDERED){

            if(rank==0)printf("Run in order mode.\n");
            ordered_update(partial_grid, k, n, s);

        }else if (e == STATIC){

            if(rank==0)printf("Run in static mode.\n");
            unsigned char* next = (unsigned char*)malloc(k*my_rows_number*sizeof(unsigned char));
            static_update(partial_grid, next, k, n, s, rank, size, my_rows_number);
            free(next);

        }else if (e == WAVE){

            if(rank==0)printf("Run in wave mode.\n");
            
            unsigned char* grid = (unsigned char*)malloc(k*k*sizeof(unsigned char));
            ///MPI_Allgatherv !!!!!
            unsigned char* next = (unsigned char*)calloc(k*k*sizeof(unsigned char), sizeof(unsigned char));
            wave_update(grid, next, k, n, s, rank, size);
            free(next);

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

