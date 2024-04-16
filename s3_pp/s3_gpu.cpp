#include <iostream>
#include <fstream>
#include <string>
#include <cuda_runtime.h>

__global__ void smith_waterman_kernel(int *dp, const char *seq1, const char *seq2, int n, int m, int match_score, int mismatch_penalty, int gap_penalty) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if (i <= n && j <= m) {
        int score_diag = dp[(i-1) * (m+1) + (j-1)] + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_penalty);
        int score_up = dp[(i-1) * (m+1) + j] + gap_penalty;
        int score_left = dp[i * (m+1) + (j-1)] + gap_penalty;
        int score = max(0, max(score_diag, max(score_up, score_left)));
        dp[i * (m+1) + j] = score;
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
    std::string seq1 = "AGTACGCA";
    std::string seq2 = "TATGC";

    int n = seq1.size();
    int m = seq2.size();
    int match_score = 2;
    int mismatch_penalty = -1;
    int gap_penalty = -1;

    int *d_dp;
    char *d_seq1, *d_seq2;

    size_t dp_size = (n + 1) * (m + 1) * sizeof(int);
    cudaMalloc(&d_dp, dp_size);
    cudaMalloc(&d_seq1, n);
    cudaMalloc(&d_seq2, m);

    cudaMemcpy(d_seq1, seq1.c_str(), n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, seq2.c_str(), m, cudaMemcpyHostToDevice);
    cudaMemset(d_dp, 0, dp_size);

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((n + threadsPerBlock.x - 1) / threadsPerBlock.x, (m + threadsPerBlock.y - 1) / threadsPerBlock.y);
    smith_waterman_kernel<<<numBlocks, threadsPerBlock>>>(d_dp, d_seq1, d_seq2, n, m, match_score, mismatch_penalty, gap_penalty);

    int *h_dp = new int[(n+1)*(m+1)];
    cudaMemcpy(h_dp, d_dp, dp_size, cudaMemcpyDeviceToHost);

    int max_score = 0;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            if (h_dp[i * (m+1) + j] > max_score) {
                max_score = h_dp[i * (m+1) + j];
            }
        }
    }

    std::cout << "Maximum score: " << max_score << std::endl;

    cudaFree(d_dp);
    cudaFree(d_seq1);
    cudaFree(d_seq2);
    delete[] h_dp;

    return 0;
}
