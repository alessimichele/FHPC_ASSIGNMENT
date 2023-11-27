void wave_update(unsigned char *grid, unsigned char* next, int k,  int n,  int s, int size, int rank);
int map_even_grid(int rand_cell_idx, int k);
unsigned int* recoverSquare(int k, int index, int radius);
void wave_update_OpenMP(unsigned char* grid, unsigned char* next, int k, int n, int s );
void wave_update_MPI(unsigned char *grid, unsigned char* next, int k, int n, int s, int size, int rank);