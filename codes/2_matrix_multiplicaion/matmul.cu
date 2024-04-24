#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N 1000
#define BLOCK_SIZE 16

__global__ void matrixMultiplyKernel(int *a, int *b, int *result, int n) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < n && col < n) {
        int sum = 0;
        for (int k = 0; k < n; k++) {
            sum += a[row * n + k] * b[k * n + col];
        }
        result[row * n + col] = sum;
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

    // Allocate memory for matrices on the host
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

    // Allocate memory for matrices on the device
    int *d_mat1, *d_mat2, *d_result;
    cudaMalloc((void **)&d_mat1, size);
    cudaMalloc((void **)&d_mat2, size);
    cudaMalloc((void **)&d_result, size);

    // Copy matrices from host to device
    cudaMemcpy(d_mat1, mat1, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_mat2, mat2, size, cudaMemcpyHostToDevice);

    // Define grid and block dimensions
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid((N + dimBlock.x - 1) / dimBlock.x, (N + dimBlock.y - 1) / dimBlock.y);

    struct timeval start, end;
    gettimeofday(&start, NULL); // Start timing

    // Launch kernel
    matrixMultiplyKernel<<<dimGrid, dimBlock>>>(d_mat1, d_mat2, d_result, N);

    // Copy result matrix from device to host
    cudaMemcpy(result, d_result, size, cudaMemcpyDeviceToHost);

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

    // Free allocated memory on the device
    cudaFree(d_mat1);
    cudaFree(d_mat2);
    cudaFree(d_result);

    // Free allocated memory on the host
    free(mat1);
    free(mat2);
    free(result);

    return 0;
}

