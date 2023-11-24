// input-output file and initialization

#include <string.h>
#include <stdlib.h>
#include <stdio.h> 

#include "io_init.h"

#define XWIDTH 256
#define YWIDTH 256
#define MAXVAL 65535


#if ((0x100 & 0xf) == 0x0)
#define I_M_LITTLE_ENDIAN 1
#define swap(mem) (( (mem) & (short int)0xff00) >> 8) +	\
  ( ((mem) & (short int)0x00ff) << 8)
#else
#define I_M_LITTLE_ENDIAN 0
#define swap(mem) (mem)
#endif



void write_pgm_image( void *grid, int maxval, int xsize, int ysize,  char *file_name)
/*
 * grid        : a pointer to the memory region that contains the grid
 * maxval       : either 255 or 65536
 * xsize, ysize : x and y dimensions of the grid
 * file_name   : the name of the file to be written
 *
 */
{
  FILE* file_stream; 

  file_stream = fopen(file_name, "w+"); 
  

  //printf("Writing file %s\n", file_name);
  //printf("file_stream address: %p\n", file_stream);

  int color_depth = 1 + ( maxval > 255 );

  if (file_stream != NULL){
    fprintf(file_stream, "P5\n# generated by\n# Alessi Michele and Carollo Marco\n%d %d\n%d\n", xsize, ysize, maxval);
    fwrite( grid, 1, xsize*ysize*color_depth, file_stream);  
    fclose(file_stream); 
  }else{
    fprintf(stderr, "Failed to open the PGM file for writing.\n");
  }
  return ;
}


void read_pgm_image( void **grid, int *maxval, int *xsize, int *ysize,  char *file_name)
/*
 * grid        : a pointer to the pointer that will contain the grid
 * maxval       : a pointer to the int that will store the maximum intensity in the grid
 * xsize, ysize : pointers to the x and y sizes
 * file_name   : the name of the file to be read
 *
 */
{
  FILE* file_stream; 
  file_stream = fopen(file_name, "r"); 

  *grid = NULL;
  *xsize = *ysize = *maxval = 0;
  
  char    MagicN[2];
  char   *line = NULL;
  size_t  k, n = 0;
  
  // get the Magic Number
  k = fscanf(file_stream, "%2s%*c", MagicN );

  // skip all the comments
  k = getline( &line, &n, file_stream);
  while ( (k > 0) && (line[0]=='#') )
    k = getline( &line, &n, file_stream);

  if (k > 0)
    {
      k = sscanf(line, "%d%*c%d%*c%d%*c", xsize, ysize, maxval);
      if ( k < 3 )
	fscanf(file_stream, "%d%*c", maxval);
    }
  else
    {
      *maxval = -1;         // this is the signal that there was an I/O error
			    // while reading the grid header
      free( line );
      return;
    }
  free( line );
  
  int color_depth = 1 + ( *maxval > 255 );
  unsigned int size = *xsize * *ysize * color_depth;
  
  if ( (*grid = (char*)malloc( size )) == NULL )
    {
      fclose(file_stream);
      *maxval = -2;         // this is the signal that memory was insufficient
      *xsize  = 0;
      *ysize  = 0;
      return;
    }
  
  if ( fread( *grid, 1, size, file_stream) != size )
    {
      free( grid );
      grid   = NULL;
      *maxval = -3;         // this is the signal that there was an i/o error
      *xsize  = 0;
      *ysize  = 0;
    }  

  fclose(file_stream);
  return;
}


void init_serial(char *file_name, int k){
  /*
    * file_name: the name of the file to be written
    * k: the size of the grid
  */

  
  unsigned char *grid;
  grid = (unsigned char*)malloc(k*k*sizeof(unsigned char));

  if (grid == NULL) {
    perror("Memory allocation error");
    return;
  }

  srand(0);

  // probability to generate a 0 or a 1
  double p = 0.2;
  // fill the grid with random values: 0 or 255
  for (int i=0; i<k*k; i++){
    grid[i] = (rand() < p * RAND_MAX) ? 0 : 255;
  }

  write_pgm_image(grid, 255, k, k, file_name);
  free(grid);

  return;
}


