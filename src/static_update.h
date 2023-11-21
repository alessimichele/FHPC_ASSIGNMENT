void static_update(unsigned char *current, unsigned char* next, int k,  int n,  int s, int rank, int size, int rows_per_process);
void static_update_OpenMP(unsigned char *current, unsigned char* next, int k,  int n,  int s);
void static_update_MPI(unsigned char *current, unsigned char* next, int k,  int n,  int s, int rank, int size, int rows_per_process);