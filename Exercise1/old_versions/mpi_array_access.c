#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// salloc -N1 -n4 -p EPYC --time=00:10:00
// [malessi@login01 scratch]$ mpicc mpi_array_access.c -o mpi_array_access
// [malessi@login01 scratch]$ mpirun mpi_array_access

int main(int argc, char** argv) {

    int N = 10;
    int rank, size;
    int *array = (int *) malloc(N * sizeof(int));

    // initialize the array
    for (int i = 0; i < N; i++) {
        array[i] = i;
    }

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    printf("Rank: %d see the array at address %p, and it is able to see all the elements\n", rank, array);
    for(int i=0; i<N; i++){
        printf("Rank: %d see the array[%d] = %d at address %p\n", rank, i, array[i], &array[i]);
    }

    
   
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0){
        printf("-------------------Scatter----------------------\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Scatter the array
    int rem = N % size;
    int sum = 0;
    int *sendcounts = (int*)malloc(size*sizeof(int));
    int *displs = (int*)malloc(size*sizeof(int));
    
    for (int i=0; i<size; i++){
        sendcounts[i] = N/size;
        if (rem > 0){
            sendcounts[i]++;
            rem--;
        }
        displs[i] = sum;
        sum += sendcounts[i];
    }


    int *recvbuf = (int *) malloc(sendcounts[rank] * sizeof(int));

    if (rank==0) {
        for (int i = 0; i < size; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatterv(array, sendcounts, displs, MPI_INT, recvbuf, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);

    // print what each process received
    printf("Rank %d has received: ", rank);
    for (int i = 0; i < sendcounts[rank]; i++) {
        printf("%d\t", recvbuf[i]);
    }
    printf("\n");

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0){
        printf("\n-------------------Do Something and Gather----------------------\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < sendcounts[rank]; i++) {
        recvbuf[i] = recvbuf[i] * rank;
    }

    // gather the results back to the root process
    MPI_Gatherv(recvbuf, sendcounts[rank], MPI_INT, array, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);


    MPI_Finalize();
    free(sendcounts);
    free(displs);
    free(recvbuf);

    // print the finale array
    if (rank == 0) {
        printf("The final array is: ");
        for (int i = 0; i < N; i++) {
            printf("%d\t", array[i]);
        }
        printf("\n");
    }
    free(array);
    return 0;
}