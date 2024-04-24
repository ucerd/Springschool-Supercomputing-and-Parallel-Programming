#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define N 1000

void matrixMultiply(int *a, int *b, int *result, int n) {
    #pragma omp parallel for collapse(1)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int sum = 0;
            for (int k = 0; k < n; k++) {
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

int main() {
    int *mat1, *mat2, *result;
    int size = N * N * sizeof(int);

    // Allocate memory for matrices
    mat1 = (int *)malloc(size);
    mat2 = (int *)malloc(size);
    result = (int *)malloc(size);
    if (mat1 == NULL || mat2 == NULL || result == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    // Initialize matrices with random values
    for (int i = 0; i < N * N; i++) {
        mat1[i] = rand() % 1000; // Random values between 0 and 999
        mat2[i] = rand() % 1000; // Random values between 0 and 999
    }

    struct timeval start, end;
    gettimeofday(&start, NULL); // Start timing

    matrixMultiply(mat1, mat2, result, N);

    gettimeofday(&end, NULL); // End timing
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;

    // Print a subset of the matrices for demonstration
    printf("Matrix 1 (subset):\n");
    printMatrix(mat1, 5, 5);
    printf("\nMatrix 2 (subset):\n");
    printMatrix(mat2, 5, 5);
    printf("\nResultant Matrix (subset):\n");
    printMatrix(result, 5, 5);
    printf("\nExecution time: %lf seconds\n", elapsed);

    // Free allocated memory
    free(mat1);
    free(mat2);
    free(result);

    return 0;
}

