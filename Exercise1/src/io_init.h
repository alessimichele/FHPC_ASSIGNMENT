#include <mpi.h>

void parallel_write(unsigned char* grid, int maxval, char* file_path, int k, int my_rows_number, int rank, int size, MPI_Comm comm);
void init_parallel(char *file_path, int k, int rank, int size, int my_rows_number);
void parallel_read( unsigned char **grid_pointer, int *maxval, int *xsize, int *ysize,  char *file_path, int my_rows_number, int rank, int size, MPI_Comm comm);