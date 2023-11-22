#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>


#include "io_init.h"
#include "ordered_update.h"

#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1

int   action = 0;
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
      fname = (char*)malloc( sizeof(optarg)+1 );
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

  if ( fname == NULL )
    {
      fname = (char*)malloc( 30 );
      sprintf(fname, "files/init/initial_gird.pgm");
    }

  //if -i is called, initialize the grid
  if ( action == INIT )
    {
      printf("initializing\n");
      // initialize the image
      init_serial(fname, k);
    }
  else if ( action == RUN )
    {
      char *grid = (char*)malloc(k*k*sizeof(char));
      void *grid_ptr = (void *)grid;
      read_pgm_image(&grid_ptr, &maxval, &k, &k, fname);
      
      printf("running\n");
      // run the simulation
      ordered_update(grid, k, n, s);
    }


 
free(fname);
  return 0;
}
