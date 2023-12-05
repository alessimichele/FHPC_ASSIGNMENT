void write_pgm_image( void *grid, int maxval, int xsize, int ysize,  char *file_name);
void read_pgm_image( void **grid, int *maxval, int *xsize, int *ysize,  char *file_name);
void init_serial( char *file_name, int k);
void parallel_write_MPI(unsigned char* grid, int maxval, char* filename, int k, int my_rows_number, MPI_Comm comm);