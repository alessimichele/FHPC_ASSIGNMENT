void wave_update(unsigned char *grid, unsigned char* next, int k,  int n,  int s, int rank, int size, int my_rows_number);
int map_even_grid(int rand_cell_idx, int k);
unsigned int* recoverSquare(int k, int index, int radius);
void wave_update_OpenMP(unsigned char* grid, unsigned char* next, int k, int n, int s, int rank, int size, int my_rows_number);
void wave_update_MPI(unsigned char *grid, unsigned char* next, int k, int n, int s, int rank, int size, int my_rows_number);