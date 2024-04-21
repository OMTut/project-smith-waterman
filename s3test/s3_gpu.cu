#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cuda_runtime.h>

__global__ void smith_waterman_kernel(int *dp, const char *seq1, const char *seq2, int n, int m, int match_score, int mismatch_penalty, int gap_penalty, int *max_score) {
    int i = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int j = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if (i > n || j > m) return;

    int index = i * (m + 1) + j;
    int score_diag = (i > 1 && j > 1) ? dp[(i - 1) * (m + 1) + (j - 1)] : ((i == 1 && j == 1) && (seq1[i-1] == seq2[j-1]) ? match_score : mismatch_penalty);
    int score_up = (i > 1) ? dp[(i - 1) * (m + 1) + j] - gap_penalty : -gap_penalty;
    int score_left = (j > 1) ? dp[i * (m + 1) + (j - 1)] - gap_penalty : -gap_penalty;

    int cell_score = max(0, max(score_diag, max(score_up, score_left)));
    int local_max_score = max(cell_score, 0);

    dp[index] = local_max_score;
    atomicMax(max_score, local_max_score);

    // Detailed diagnostics
    if (threadIdx.x == 0 || threadIdx.x == blockDim.x - 1 || threadIdx.y == 0 || threadIdx.y == blockDim.y - 1) {
        printf("Boundary Thread (%d,%d) -> idx=%d, diag=%d, up=%d, left=%d, cell=%d, local_max=%d, global_max=%d\n",
               i, j, index, score_diag, score_up, score_left, local_max_score, *max_score);
    }
}



std::string read_fasta_block(std::ifstream& file, size_t block_size) {
    std::string block, line;
    while (block.size() < block_size && std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;  // Skip headers or empty lines
        block += line;
    }
    return block;
}

int smith_waterman_cuda(const std::string& seq1, const std::string& seq2, int match_score, int mismatch_penalty, int gap_penalty) {
    int n = seq1.length();
    int m = seq2.length();
    int *dp, *max_score;
    char *d_seq1, *d_seq2;

    cudaMalloc((void**)&dp, (n+1) * (m+1) * sizeof(int));
    cudaMalloc((void**)&max_score, sizeof(int));
    cudaMemset(max_score, 0, sizeof(int));  // Ensure max_score is zeroed before use

    cudaMalloc((void**)&d_seq1, (n+1) * sizeof(char));
    cudaMalloc((void**)&d_seq2, (m+1) * sizeof(char));
    cudaMemcpy(d_seq1, seq1.c_str(), n+1, cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, seq2.c_str(), m+1, cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((m + threadsPerBlock.x - 1) / threadsPerBlock.x, (n + threadsPerBlock.y - 1) / threadsPerBlock.y);

    smith_waterman_kernel<<<numBlocks, threadsPerBlock>>>(dp, d_seq1, d_seq2, n, m, match_score, mismatch_penalty, gap_penalty, max_score);
    cudaDeviceSynchronize();  // Ensure all operations are completed

    int result;
    cudaMemcpy(&result, max_score, sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(dp);
    cudaFree(d_seq1);
    cudaFree(d_seq2);
    cudaFree(max_score);

    return result;
}

int main(int argc, char* argv[]) {
    std::string path1 = argv[1];
    std::string path2 = argv[2]; 
    std::ifstream file1(path1), file2(path2);
    const size_t block_size = 10000;  // Adjust block size based on your memory constraints

    if (!file1.is_open() || !file2.is_open()) {
        std::cerr << "Error opening files!" << std::endl;
        return -1;
    }

    std::vector<std::string> seq2_blocks;
    std::string seq2_block;
    while (!(seq2_block = read_fasta_block(file2, block_size)).empty()) {
        seq2_blocks.push_back(seq2_block);
    }
    file2.close();

    std::string seq1_block;
    while (!(seq1_block = read_fasta_block(file1, block_size)).empty()) {
        for (const auto& block : seq2_blocks) {
            int score = smith_waterman_cuda(seq1_block, block, 2, -1, -1);
            std::cout << "Alignment score: " << score << std::endl;
        }
    }

    file1.close();
    return 0;
}
