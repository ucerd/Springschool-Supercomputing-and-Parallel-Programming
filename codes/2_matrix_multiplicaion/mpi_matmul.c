#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define N 1000 // Debugging with smaller matrix size

void matrixMultiply(int *a, int *b, int *result, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int sum = 0;
            for (int k = 0; k < n; k++) {
               // printf("i=%d, j=%d, k=%d\n", i, j, k); // Debug print statement
                sum += a[i * n + k] * b[k * n + j];
            }
            result[i * n + j] = sum;
        }
    }
}

void printMatrix(int *mat, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%d\t", mat[i * cols + j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    int *mat1, *mat2, *result;
    int size = N * N * sizeof(int);
    int rank, numprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    // Allocate memory for matrices on each process
    mat1 = (int *)malloc(size);
    mat2 = (int *)malloc(size);
    result = (int *)malloc(size);
    if (mat1 == NULL || mat2 == NULL || result == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Initialize matrices with random values on rank 0 process
    if (rank == 0) {
        for (int i = 0; i < N * N; i++) {
            mat1[i] = rand() % 1000; // Random values between 0 and 999
            mat2[i] = rand() % 1000; // Random values between 0 and 999
        }
    }

    MPI_Bcast(mat1, N * N, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast mat1 to all processes
    MPI_Bcast(mat2, N * N, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast mat2 to all processes

    struct timeval start, end;
    gettimeofday(&start, NULL); // Start timing

    matrixMultiply(mat1, mat2, result, N);

    gettimeofday(&end, NULL); // End timing
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;

    // Gather results from all processes to rank 0
    int *finalResult = NULL;
    if (rank == 0) {
        finalResult = (int *)malloc(N * N * sizeof(int));
        if (finalResult == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    MPI_Gather(result, N * N / numprocs, MPI_INT, finalResult, N * N / numprocs, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Print a subset of the matrices for demonstration
        printf("\nExecution time: %lf seconds\n", elapsed);

        free(finalResult);
    }

    // Free allocated memory
    free(mat1);
    free(mat2);
    free(result);

    MPI_Finalize();

    return 0;
}

