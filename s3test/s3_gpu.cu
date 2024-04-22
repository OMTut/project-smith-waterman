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

// Kernel to initialize the DP table
__global__ void init_dp(int *d_dp, int n, int m) {
    int row = blockIdx.x * blockDim.x + threadIdx.x; // Row index
    int col = blockIdx.y * blockDim.y + threadIdx.y; // Column index

    if (row <= n && col <= m) {
        d_dp[row * (m + 1) + col] = 0; // Initialize to zero
        printf("GPU: %d, Row: %d, Col: %d, Initial DP Value: %d\n", blockIdx.z, row, col, d_dp[row * (m + 1) + col]); // Print initial DP value
    }
}

// Kernel for the Smith-Waterman algorithm using diagonal stripes
__global__ void smith_waterman(int *d_dp, const char *d_seq1, const char *d_seq2, int n, int m) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    printf("Seq1: %c, Seq2: %c\n", d_seq1[row - 1], d_seq2[col - 1]);

    for (int sum = 2; sum <= n + m; sum++) {
        if (row > 0 && row <= n && col > 0 && col <= m && row + col == sum) {
            int diag_score = d_dp[(row - 1) * (m + 1) + (col - 1)];
            int top_score = d_dp[(row - 1) * (m + 1) + col];
            int left_score = d_dp[row * (m + 1) + (col - 1)];

            int match = diag_score + ((d_seq1[row - 1] == d_seq2[col - 1]) ? MATCH_SCORE : -MISMATCH_PENALTY);
            int delete_score = top_score - GAP_PENALTY;
            int insert_score = left_score - GAP_PENALTY;

            int final_score = max(max(match, delete_score), max(insert_score, 0));

            d_dp[row * (m + 1) + col] = final_score;
            printf("Match: %d, Delete: %d, Insert: %d\n", match, delete_score, insert_score);
            printf("Row: %d, Col: %d, Score: %d\n", row, col, final_score);
        }
        __syncthreads();
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

    //count gpus and divide the sequence into parts
    int num_gpus = getNumberofGPUS();
    std::cout << "Found " << num_gpus << " CUDA device(s)" << std::endl;
    int part_size = n / num_gpus;

    int *dp = new int[(n + 1) * (m + 1)](); // Initialize DP table in host memory

    for (int gpu = 0; gpu < num_gpus; gpu++) {
        cudaSetDevice(gpu);
        int start = gpu * part_size;
        int end = start + part_size;
        int size = end - start;
        std::cout << "GPU " << gpu << " will process sequence from " << start << " to " << end << std::endl;

        // Initialize the dp matrix with zeroes
        int *dp = new int[(size + 1) * (m + 1)](); // Initialize with zeros
        int *d_dp;
        cudaMalloc(&d_dp, (size + 1) * (m + 1) * sizeof(int));
        cudaMemset(d_dp, 0, (size + 1) * (m + 1) * sizeof(int));

        char *d_seq1, *d_seq2;
        cudaMalloc(&d_seq1, size * sizeof(char));
        cudaMalloc(&d_seq2, m * sizeof(char));
        cudaMemcpy(d_seq1, seq1.substr(start, size).c_str(), size * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(d_seq2, seq2.c_str(), m * sizeof(char), cudaMemcpyHostToDevice);

        dim3 threadsPerBlock(16, 16); // Adjusted block size for computation
        dim3 numBlocks((size + 15) / 16, (m + 15) / 16); // Grid size for computation
        smith_waterman<<<numBlocks, threadsPerBlock>>>(d_dp, d_seq1, d_seq2, size, m);

        cudaDeviceSynchronize(); // Ensure computation is complete

        int *temp_dp = new int[(size + 1) * (m + 1)];
        cudaMemcpy(temp_dp, d_dp, (size + 1) * (m + 1) * sizeof(int), cudaMemcpyDeviceToHost); // Copy DP table back
        // Accumulate DP table across all GPUs
        for (int i = 0; i <= size; i++) {
            for (int j = 0; j <= m; j++) {
                dp[i * (m + 1) + j] += temp_dp[i * (m + 1) + j];
            }
        }

        delete[] temp_dp;
        cudaFree(d_dp);
        cudaFree(d_seq1);
        cudaFree(d_seq2);
    }

    //merge results
    // Find the maximum score and its position in the entire DP matrix
    int max_score = 0, max_i = 0, max_j = 0;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            int score = dp[i * (m + 1) + j];
            if (score > max_score) {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Perform the traceback only once, starting at the cell with the maximum score in the entire DP matrix
    int alignment = traceback(dp, seq1.c_str(), seq2.c_str(), max_i, max_j);

    std::cout << "Alignment: " << alignment << std::endl;

    delete[] dp;

    return 0;
}

