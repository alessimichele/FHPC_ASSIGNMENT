#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>



void parallel_write_MPI(unsigned char* grid, int maxval, char* filename, int k, int my_rows_number, MPI_Comm comm) {
    
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Status status;
    //initialize the variable header_offset in each process
    MPI_Offset header_offset;
    if (rank==0){
        FILE* file_stream; 

        file_stream = fopen(filename, "w+"); //file is going to be either created or overwritten,
        // and the file pointer is positioned at the beginning of the file,
        //so the function ftell will return 0, which is the initial position of the file pointer

        int color_depth = 1 + ( maxval > 255 );
        
        if (file_stream != NULL){
            fprintf(file_stream, "P5\n# generated by\n# Alessi Michele and Carollo Marco\n%d %d\n%d\n", k, k, maxval);
            header_offset = ftell(file_stream);
            fclose(file_stream); 
        }else{
            fprintf(stderr, "Failed to open the PGM file for writing.\n");
        }
    }
    MPI_Barrier(comm);  
    
    MPI_File file;
    MPI_File_open(comm, filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

    // Now broadcast the initial position (from rank 0, that got it with ftell) to all other processes
    MPI_Bcast(&header_offset, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    MPI_Barrier(comm);
    //let's calculate the offset, first we need an array with the number of rows for each process
    int* rows_per_process = (int*)malloc(size*sizeof(int));
    for (int i=0; i<size; i++){
        rows_per_process[i] = (i<(k%size)) ? k/size +1 : k/size;
    }
    int* offset_arr = (int*)malloc(size*sizeof(int));
    offset_arr[0] = 0;
    for (int i=1; i<size; i++){
        offset_arr[i] = offset_arr[i-1] + rows_per_process[i-1]*k;
    }
    free(rows_per_process);

    MPI_Offset offset = offset_arr[rank]*sizeof(unsigned char)+header_offset;

    MPI_File_seek(file, offset, MPI_SEEK_SET);
    MPI_Barrier(comm);
    MPI_File_write(file, grid, my_rows_number * k, MPI_UNSIGNED_CHAR, &status);

    MPI_File_close(&file);
    free(offset_arr);
    return;
}


void init_parallel(char *file_name, int k, int rank, int size){

    int my_rows_number = (rank<(k%size)) ? k/size+1 : k/size;
    unsigned char *grid = (unsigned char*)malloc(k*my_rows_number*sizeof(unsigned char));
    
    if (grid == NULL) {
      perror("Memory allocation error");
      return;
    }

    unsigned int seed = clock();
    seed = seed * rank;
   
    for (int i=0; i<my_rows_number*k; i++){
        grid[i] = rand_r(&seed)%2 * 255;
    }
    ////print the grid
    //for (int i=0; i<my_rows_number; i++){
    //    for (int j=0; j<k; j++){
    //        printf("%d ", grid[i*k+j]);
    //    }
    //    printf("\n");
    //}

    parallel_write_MPI(grid, 255, file_name, k, my_rows_number, MPI_COMM_WORLD);
    free(grid);
    return;
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


void parallel_read_MPI( unsigned char **grid_pointer, int *maxval, int *xsize, int *ysize,  char *file_name, int my_rows_number, MPI_Comm comm) {
    
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Status status;
    //initialize the variable header_offset in each process
    MPI_Offset header_offset;
    if (rank==0){
       
        FILE* file_stream; 
        file_stream = fopen(file_name, "r"); 
        
        *xsize = *ysize = *maxval = 0;
        
        char    MagicN[2];
        char   *line = NULL;
        size_t  m, n = 0; // m is the number of characters read, n is the size of the buffer
        
        // get the Magic Number
        m = fscanf(file_stream, "%2s%*c", MagicN );
        header_offset+=m+1; //+1 because of the \n
        // skip all the comments
        m = getline( &line, &n, file_stream);
        header_offset+=m;
        while ( (m > 0) && (line[0]=='#') )
            m = getline( &line, &n, file_stream);
            header_offset+=m;   
        if (m > 0)
          {
            m = sscanf(line, "%d%*c%d%*c%d%*c", xsize, ysize, maxval);
            //maybe fix HEADER_OFFSET
            if ( m < 3 )
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
        
    }
    MPI_Barrier(comm);  
    
    MPI_File file;
    MPI_File_open(comm, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    // Now broadcast the initial position (from rank 0, that got it with ftell) to all other processes
    MPI_Bcast(&header_offset, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    //broadcast also xsize
    MPI_Bcast(xsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    *grid_pointer = (unsigned char*)malloc(*xsize*my_rows_number*sizeof(unsigned char));
    MPI_Barrier(comm);
    //let's calculate the offset, first we need an array with the number of rows for each process
    int* rows_per_process = (int*)malloc(size*sizeof(int));
    for (int i=0; i<size; i++){
        rows_per_process[i] = (i<(*xsize%size)) ? *xsize/size +1 : *xsize/size;
    }
    int* offset_arr = (int*)malloc(size*sizeof(int));
    offset_arr[0] = 0;
    for (int i=1; i<size; i++){
        offset_arr[i] = offset_arr[i-1] + rows_per_process[i-1]*(*xsize);
    }
    free(rows_per_process);

    MPI_Offset offset = offset_arr[rank]*sizeof(unsigned char)+header_offset;

    MPI_File_seek(file, offset, MPI_SEEK_SET);
    MPI_Barrier(comm);
    MPI_File_read(file, *grid_pointer, my_rows_number * (*xsize), MPI_UNSIGNED_CHAR, &status);


    MPI_File_close(&file);
    free(offset_arr);

    return;

}