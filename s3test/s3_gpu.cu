#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cuda_runtime.h>

#define MATCH_SCORE     2
#define MISMATCH_PENALTY 1
#define GAP_PENALTY    1

int getNumberofGPUS() {
    int num_gpus = 0;
    cudaGetDeviceCount(&num_gpus);
    return num_gpus;
}

__global__ void init_dp(int *d_dp, int n, int m) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    int col = blockIdx.y * blockDim.y + threadIdx.y;

    if (row <= n && col <= m) {
        d_dp[row * (m + 1) + col] = 0;
    }
}

__global__ void smith_waterman(int *d_dp, const char *d_seq1, const char *d_seq2, int n, int m) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row > 0 && row <= n && col > 0 && col <= m) {
        int diag_score = d_dp[(row - 1) * (m + 1) + (col - 1)];
        int top_score = d_dp[(row - 1) * (m + 1) + col];
        int left_score = d_dp[row * (m + 1) + (col - 1)];

        int match = diag_score + ((d_seq1[row - 1] == d_seq2[col - 1]) ? MATCH_SCORE : -MISMATCH_PENALTY);
        int delete_score = top_score - GAP_PENALTY;
        int insert_score = left_score - GAP_PENALTY;

        int final_score = max(max(match, delete_score), max(insert_score, 0));

        d_dp[row * (m + 1) + col] = final_score;
    }
}

std::string read_sequence(const std::string &filename) {
    std::ifstream file(filename);
    std::string line, sequence;
    bool first_line = true;
    while (std::getline(file, line)) {
        if (first_line) {
            first_line = false;
            continue;
        }
        sequence += line;
    }
    return sequence;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <seq1.fasta> <seq2.fasta>" << std::endl;
        return -1;
    }

    int num_gpus = getNumberofGPUS();

    cudaStream_t *streams = new cudaStream_t[num_gpus];
    for (int i = 0; i < num_gpus; i++) {
        cudaStreamCreate(&streams[i]);
    }

    std::cout << "Running Smith-Waterman algorithm on GPU" << std::endl;
    std::cout << "Getting sequences from FASTA files" << std::endl;
    std::string seq1 = read_sequence(argv[1]);
    std::string seq2 = read_sequence(argv[2]);

    double available_memory = 0.8 * 12 * 1024 * 1024 * 1024;
    double max_chunk_size = sqrt(available_memory / 4);
    size_t chunk_size = min(static_cast<size_t>(max_chunk_size), max(seq1.size(), seq2.size()));
    std::cout << "Chunk size: " << chunk_size <<std::endl;
    int num_chunks_processed = 0;
    int max_score = 0;
    for (int i = 0; i < num_gpus; i++) {
        cudaSetDevice(i);

        char* d_seq1;
        char* d_seq2;
        int* d_dp;
        cudaMalloc(&d_seq1, chunk_size * sizeof(char));
        cudaMalloc(&d_seq2, chunk_size * sizeof(char));
        cudaMalloc(&d_dp, chunk_size * chunk_size * sizeof(int));

        cudaMemcpy(d_seq1, seq1.c_str(), chunk_size * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(d_seq2, seq2.c_str(), chunk_size * sizeof(char), cudaMemcpyHostToDevice);

        dim3 threadsPerBlock(16, 16);
        dim3 numBlocks((chunk_size + threadsPerBlock.x - 1) / threadsPerBlock.x, (chunk_size + threadsPerBlock.y - 1) / threadsPerBlock.y);

        init_dp<<<numBlocks, threadsPerBlock, 0, streams[i]>>>(d_dp, chunk_size, chunk_size);
        smith_waterman<<<numBlocks, threadsPerBlock, 0, streams[i]>>>(d_dp, d_seq1, d_seq2, chunk_size, chunk_size);

        cudaDeviceSynchronize();

        int* scores = new int[chunk_size * chunk_size];
        cudaMemcpy(scores, d_dp, chunk_size * chunk_size * sizeof(int), cudaMemcpyDeviceToHost);

        // Find the maximum score
        max_score = scores[0];
        for (int j = 0; j < chunk_size * chunk_size; j++) {
            if (scores[j] > max_score) {
                max_score = scores[j];
            }
        }
        num_chunks_processed++;

        delete[] scores;
        cudaFree(d_dp);
        cudaFree(d_seq1);
        cudaFree(d_seq2);
    }
    std::cout << "Maximum alignment score: " << max_score << std::endl;
    std::cout << "Number of chunks processed: " << num_chunks_processed << std::endl;
    for (int i = 0; i < num_gpus; i++) {
        cudaStreamDestroy(streams[i]);
    }
    delete[] streams;

    return 0;
}