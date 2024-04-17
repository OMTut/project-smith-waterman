#include <iostream>
#include <fstream>
#include <string>
#include <cuda_runtime.h>

__global__ void smith_waterman_kernel(int *dp, const char *seq1, const char *seq2, int n, int m, int match_score, int mismatch_penalty, int gap_penalty) {
    int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if (i <= n && j <= m) {
        int index = i * (m + 1) + j;
        int score_diag = dp[(i-1) * (m+1) + (j-1)];
        int score_up = dp[(i-1) * (m+1) + j];
        int score_left = dp[i * (m+1) + (j-1)];

        if (seq1[i-1] == seq2[j-1]) {
            score_diag += match_score;
        } else {
            score_diag -= mismatch_penalty;
        }
        score_up -= gap_penalty;
        score_left -= gap_penalty;

        int score = max(0, max(score_diag, max(score_up, score_left)));
        dp[index] = score;
    }
}


std::string read_fasta_block(std::ifstream& file, size_t block_size) {
    std::string block, line;
    while (block.size() < block_size && std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        block += line;
    }
    return block;
}

int main() {
   // the max score with this is 8
    std::string seq1 = "GATC";
    std::string seq2 = "GATC";

    int n = seq1.size();
    int m = seq2.size();
    int match_score = 2;
    int mismatch_penalty = 1;
    int gap_penalty = 1;

    int *d_dp;
    char *d_seq1, *d_seq2;

    size_t dp_size = (n + 1) * (m + 1) * sizeof(int);
    cudaMalloc(&d_dp, dp_size);
    cudaMalloc(&d_seq1, n + 1);
    cudaMalloc(&d_seq2, m + 1);

    cudaMemcpy(d_seq1, seq1.c_str(), n + 1, cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, seq2.c_str(), m + 1, cudaMemcpyHostToDevice);
    cudaMemset(d_dp, 0, dp_size);

   // Execute the kernel
   // Assume maximum sizes you would like per block
   dim3 threadsPerBlock(16, 16);  // 16x16 threads per block is common

   // Calculate number of blocks needed in each dimension
   dim3 numBlocks((m + threadsPerBlock.x - 1) / threadsPerBlock.x, 
                  (n + threadsPerBlock.y - 1) / threadsPerBlock.y);

   // Kernel launch
   smith_waterman_kernel<<<numBlocks, threadsPerBlock>>>(d_dp, d_seq1, d_seq2, n, m, match_score, mismatch_penalty, gap_penalty);


   cudaError_t error = cudaGetLastError();
   if (error != cudaSuccess) {
      std::cerr << "CUDA error: " << cudaGetErrorString(error) << std::endl;
      return -1;  // Or handle the error appropriately
   }
   cudaDeviceSynchronize();

    int *h_dp = new int[(n+1)*(m+1)];
    cudaMemcpy(h_dp, d_dp, dp_size, cudaMemcpyDeviceToHost);

   //  int max_score = 0;
   //  for (int i = 0; i <= n; ++i) {
   //      for (int j = 0; j <= m; ++j) {
   //          if (h_dp[i * (m+1) + j] > max_score) {
   //              max_score = h_dp[i * (m+1) + j];
   //          }
   //      }
   //  }

   int max_score = 0;
   for (int i = 0; i <= n; i++) {
      for (int j = 0; j <= m; j++) {
         max_score = max(max_score, h_dp[i * (m+1) + j]);
      }
   }

    std::cout << "Maximum score: " << max_score << std::endl;

    // Clean up memory allocations
    cudaFree(d_dp);
    cudaFree(d_seq1);
    cudaFree(d_seq2);
    delete[] h_dp;

    return 0;
}
