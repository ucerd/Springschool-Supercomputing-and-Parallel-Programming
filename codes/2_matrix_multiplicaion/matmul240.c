#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 480 

void printMatrix(int mat[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d\t", mat[i][j]);
        }
        printf("\n");
    }
}

void matrixMultiply(int mat1[N][N], int mat2[N][N], int result[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i][j] = 0;
            for (int k = 0; k < N; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

int main(int argc, char** argv) {
    int rank, size;
    int mat1[N][N];
    int mat2[N][N];
    int result[N][N];

    // Initialize matrices mat1 and mat2
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            mat1[i][j] = i * N + j + 1; // Fill with consecutive numbers
            mat2[i][j] = (N * N) - (i * N + j); // Fill with numbers in reverse order
        }
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 240) {
        printf("Error: This program requires exactly 240 processes.\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Start the timer
    double startTime = MPI_Wtime();

    matrixMultiply(mat1, mat2, result);

    // End the timer
    double endTime = MPI_Wtime();

    if (rank == 0) {
        printf("Matrix 1:\n");
        printMatrix(mat1);
        printf("\nMatrix 2:\n");
        printMatrix(mat2);
        printf("\nResultant Matrix:\n");
        printMatrix(result);
        printf("\nTime elapsed: %lf seconds\n", endTime - startTime);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}

