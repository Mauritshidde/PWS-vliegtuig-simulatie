#include <iostream>

__global__ void addvec(int N, int M, int B, double *array, double *sum) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    // int j = blockIdx.y * blockDim.y + threadIdx.y;
    // int index = j + i * M;
    if (i < N * M * B) {
        // sum[0] += array[i];
        __syncthreads();
        array[i] = 3;
        // printf("%d ",i);
    } else {
        __syncthreads();
        sum[0] = i;
    }
    
}

int main() {
    int N = 40;
    int M = 40;
    int B = 40;
    int block_size = 256;
    int grid_size = ((N * M * B + 255) / block_size);
    // dim3 threadsPerBlock((N+255)/256, (M+255)/256, (B+255)/256);

    double *array = (double*)malloc(N * M * B * sizeof(double));
    double *sum = (double*)malloc(sizeof(double));
    sum[0] = 2;

    for (int i=0; i < N; i++) {
        for (int j=0; j < M; j++) {
            for (int k=0; k < B; k++) {
                array[k + i*M + j * B] = 1;
            }
        }
    } 

    double *array_p;
    double *sum_p;

    cudaMalloc(&array_p, N * M * B *sizeof(double));
    cudaMemcpy(array_p, array, N * M * B * sizeof(double), cudaMemcpyHostToDevice);
    
    cudaMalloc(&sum_p, sizeof(double));
    cudaMemcpy(sum_p, sum, sizeof(double), cudaMemcpyHostToDevice);

    addvec<<<grid_size, block_size>>>(N, M, B, array_p, sum_p);

    cudaMemcpy(array, array_p, N * M * B * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(sum, sum_p, sizeof(double), cudaMemcpyDeviceToHost);


    // for (int i=0; i < N * M; i++) {
    //     std::cout << array[i] << " ";
    // }
    for (int i=0; i < N; i++) {
        for (int j=0; j < M; j++) {
            for (int k=0; k < B; k++) {
                std::cout << array[j + i*M] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    } 
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << sum[0] << " the sum " << std::endl;

    cudaFree(array_p);
    free(array);
    
    cudaFree(sum_p);
    free(sum);

    return 0;
}