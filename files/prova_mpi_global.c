#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>


void main(int argc, char **argv)
{  
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("Hello world from process %d of %d\n", rank, size);
    
    #pragma omp parallel
    {  
        #pragma omp master
        {
            int nthreads = omp_get_num_threads();
            printf("There are %d threads\n", nthreads);
        }
        int my_thread_id = omp_get_thread_num();

        printf( "\tgreetings from thread num %d of %d, in process %d\n", my_thread_id, omp_get_num_threads(), rank);
    }

    
    MPI_Finalize();
}
