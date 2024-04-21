#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cuda_runtime.h>

#define MATCH_SCORE     2
#define MISMATCH_PENALTY 1
#define GAP_PENALTY    1

// Kernel to initialize the DP table
__global__ void init_dp(int *d_dp, int n, int m) {
    int row = blockIdx.x * blockDim.x + threadIdx.x; // Row index
    int col = blockIdx.y * blockDim.y + threadIdx.y; // Column index

    if (row <= n && col <= m) {
        d_dp[row * (m + 1) + col] = 0; // Initialize to zero
    }
}

// Kernel for the Smith-Waterman algorithm
__global__ void smith_waterman(int *d_dp, const char *d_seq1, const char *d_seq2, int n, int m, int start, int sum) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x; // Index along the diagonal
    if (idx >= sum - start + 1) return; // Out of bounds

    int row = start + idx;
    int col = sum - row;

    if (row > 0 && row <= n && col > 0 && col <= m) {
        int diag_score = d_dp[(row - 1) * (m + 1) + (col - 1)];
        int top_score = d_dp[(row - 1) * (m + 1) + col];
        int left_score = d_dp[row * (m + 1) + (col - 1)];

        printf("Thread [%d, %d] read scores: diag=%d, top=%d, left=%d\n", threadIdx.x, threadIdx.y, diag_score, top_score, left_score);

        int match = diag_score + ((d_seq1[row - 1] == d_seq2[col - 1]) ? MATCH_SCORE : -MISMATCH_PENALTY);
        int delete_score = top_score - GAP_PENALTY;
        int insert_score = left_score - GAP_PENALTY;

        printf("Thread [%d, %d] calculated scores: match=%d, delete=%d, insert=%d\n", threadIdx.x, threadIdx.y, match, delete_score, insert_score);

        int final_score = max(max(match, delete_score), max(insert_score, 0));

        printf("Thread [%d, %d] calculated final score: %d\n", threadIdx.x, threadIdx.y, final_score);

        d_dp[row * (m + 1) + col] = final_score;
    }
    __syncthreads(); // Ensure all threads have updated the DP table before printing

    if (threadIdx.x == 0 && threadIdx.y == 0) { // Only one thread prints the DP table
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                printf("%d ", d_dp[i * (m + 1) + j]);
            }
            printf("\n");
        }
    }
}

// Read sequence from FASTA file
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

// Traceback function to build alignment
int traceback(int *dp, const std::string &seq1, const std::string &seq2, int n, int m) {
    std::string alignment1, alignment2;
    int max_score = 0, max_i = 1, max_j = 1;

    // Find the position of the maximum score
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            if (dp[i * (m + 1) + j] > max_score) {
                max_score = dp[i * (m + 1) + j];
                max_i = i;
                max_j = j;
            }
        }
    }
    printf("Starting traceback at cell (%d, %d) with score %d\n", max_i, max_j, max_score);
    int i = max_i;
    int j = max_j;

    // Traceback to build the alignment
    while (i > 0 && j > 0) {
        int current_score = dp[i * (m + 1) + j];

        int diag_score = dp[(i - 1) * (m + 1) + (j - 1)];
        int up_score = dp[(i - 1) * (m + 1) + j];
        int left_score = dp[i * (m + 1) + (j - 1)];

        printf("Calculating cell (%d, %d) with score %d\n", i, j, current_score);

        if (current_score == diag_score + ((seq1[i - 1] == seq2[j - 1]) ? MATCH_SCORE : -MISMATCH_PENALTY)) {
            alignment1 += seq1[i - 1];
            alignment2 += seq2[j - 1];
            i--;
            j--;
        } else if (current_score == up_score - GAP_PENALTY) {
            alignment1 += seq1[i - 1];
            alignment2 += '-';
            i--;
        } else if (current_score == left_score - GAP_PENALTY) {
            alignment1 += '-';
            alignment2 += seq2[j - 1];
            j--;
        } else {
            break; // Safety to avoid infinite loops
        }
    }

    std::reverse(alignment1.begin(), alignment1.end());
    std::reverse(alignment2.begin(), alignment2.end());
    return max_score;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <seq1.fasta> <seq2.fasta>" << std::endl;
        return -1;
    }

    std::string seq1 = read_sequence(argv[1]);
    std::string seq2 = read_sequence(argv[2]);
    int n = seq1.size();
    int m = seq2.size();
    //std::cout << "Sequence 1: " << seq1 << " - Sequence 1 size: " << n << std::endl;
    //std::cout << "Sequence 2: " << seq2 << " - Sequence 2 size: " << m << std::endl;

    int *dp = new int[(n + 1) * (m + 1)](); // Initialize with zeros
    int *d_dp;
    cudaMalloc(&d_dp, (n + 1) * (m + 1) * sizeof(int));

    char *d_seq1, *d_seq2;
    cudaMalloc(&d_seq1, n * sizeof(char));
    cudaMalloc(&d_seq2, m * sizeof(char));
    cudaMemcpy(d_seq1, seq1.c_str(), n * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq2, seq2.c_str(), m * sizeof(char), cudaMemcpyHostToDevice);

    dim3 threadsPerBlockInit(16, 16); // Adjusted block size for initialization
    dim3 numBlocksInit((n + 15) / 16, (m + 15) / 16); // Grid size for initialization
    //printf("Launching init_dp kernel with %d blocks of %d threads\n", numBlocksInit.x * numBlocksInit.y, threadsPerBlockInit.x * threadsPerBlockInit.y);
    init_dp<<<numBlocksInit, threadsPerBlockInit>>>(d_dp, n, m);

    cudaDeviceSynchronize(); // Wait for the GPU to finish
    //printf("Finished init_dp kernel\n");

    dim3 threadsPerBlock(16, 16); // Adjusted block size for computation
    for (int sum = 0; sum <= n + m; sum++) {
        int start = max(0, sum - m);
        int end = min(sum, n);
        int numDiagonalCells = end - start + 1;
        int numBlocks = (numDiagonalCells + 16 - 1) / 16;
        smith_waterman<<<numBlocks, threadsPerBlock>>>(d_dp, d_seq1, d_seq2, n, m, start, sum);
        cudaDeviceSynchronize(); // Wait for all threads to finish
    }

    cudaDeviceSynchronize(); // Ensure computation is complete

    cudaMemcpy(dp, d_dp, (n + 1) * (m + 1) * sizeof(int), cudaMemcpyDeviceToHost); // Copy DP table back

    int alignment = traceback(dp, seq1, seq2, n, m);
    std::cout << "Alignment: " << alignment << std::endl;

    delete[] dp;
    cudaFree(d_dp);
    cudaFree(d_seq1);
    cudaFree(d_seq2);

    return 0;
}
