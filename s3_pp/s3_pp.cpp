#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <omp.h>  // Include OpenMP header

std::string read_fasta_block(std::ifstream& file, size_t block_size) {
    std::string block, line;
    while (block.size() < block_size && std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;  // Skip headers or empty lines
        block += line;
    }
    return block;
}

int smith_waterman(const std::string& seq1, const std::string& seq2, int match_score, int mismatch_penalty, int gap_penalty) {
    int n = seq1.length();
    int m = seq2.length();
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    int max_score = 0;

    /*************************************
     * Parallelized the inner loop. WHen I try to do the out loop,
     * the results are incoorect. I think it is because the max_score
     * is not being updated correctly or a race condition.
    */
    for (int i = 1; i <= n; ++i) {
        #pragma omp parallel for reduction(max:max_score)
        for (int j = 1; j <= m; ++j) {
            int score_diagonal = dp[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_penalty);
            int score_up = dp[i-1][j] + gap_penalty;
            int score_left = dp[i][j-1] + gap_penalty;
            dp[i][j] = std::max(0, std::max(score_diagonal, std::max(score_up, score_left)));
            if (dp[i][j] > max_score) max_score = dp[i][j]; // Use reduction to handle max_score safely
        }
    }
    return max_score;
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
            int score = smith_waterman(seq1_block, block, 2, -1, -1);
            std::cout << "Alignment score: " << score << std::endl;
        }
    }

    file1.close();
    return 0;
}
