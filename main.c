#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>


#include "io_init.h"
#include "ordered_update.h"
#include "ordered_update_finite.h"
#include "static_update.h"
#include "wave_update.h"

#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1
#define ORDERED_FINITE 2

int   action = INIT;
int   k      = K_DFLT;
int   e      = ORDERED;
int   n      = 10000;
int   s      = 1;
char *fname  = NULL;
int maxval = 65535;

int main ( int argc, char **argv )
{
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
      fname = (char*)malloc(strlen(optarg)*sizeof(char) + 1);
      if (fname == NULL) {
        perror("Memory allocation error");
        return 1;
      }
      strcpy(fname, optarg);
      break;

    case 'n':
      n = atoi(optarg); break;

    case 's':
      s = atoi(optarg); break;

    default :
      printf("argument -%c not known\n", c ); break;
    }
  }
  
  // Where the initial grid is stored
  char *path = (char*)malloc(11*sizeof(char) + 1);
  strcpy(path, "files/init/");
  char *file_path = (char*)malloc((strlen(path) + strlen(fname))*sizeof(char) + 1);
  strcpy(file_path, path);
  
  if (fname != NULL) {
    strcat(file_path, fname);
} else {
    perror("Filename is not provided. Please provide a filename with -f option. This will be the name of the file containing the initial grid.\n");
    return 1; // return with error code
}



  // if -i is called, initialize the grid
  if ( action == INIT )
    { 
      printf("Initializing...\n");
      init_serial(file_path, k);
    }else{ 
    
    if (action != RUN) {
      perror("Please provide an action. Either -i for initialization or -r for running the simulation.\n");
      return 1; 
    }

    
    unsigned char* grid;
    
    read_pgm_image((void**)&grid, &maxval, &k, &k, file_path);
    printf("The initial grid has been read.\n");
   

  if ( e == 0 )
    {
      printf("Running in ordered mode...\n");

      ordered_update(grid, k, n, s);

      printf("Completed.\n");
    }else if(e == 1){
      printf("Running in static mode...\n");
      
      unsigned char* next = (unsigned char*)malloc(k*k*sizeof(unsigned char));
      static_update_OpenMP(grid, next, k, n, s);
      free(next);

      printf("Completed.\n");
    }else if(e == 2){
      printf("Running in ordered finite mode...\n");

      ordered_update_finite(grid, k, n, s);

      printf("Completed.\n");
    }else{
      printf("Running in wave mode...\n");
      
      unsigned char* next = (unsigned char*)malloc(k*k*sizeof(unsigned char));
      wave_update(grid, next, k, n, s);
      free(next);
      
      printf("Completed.\n");
    }
    free(grid);
}
  free(fname);
  free(path);
  free(file_path);
  return 0;
}
