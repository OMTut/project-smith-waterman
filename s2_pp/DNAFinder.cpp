#include "DNAFinder.h"
#include <iostream>
#include <fstream>
#include <string>
#include <future>
#include <vector>
#include <atomic>

using std::cout, std::cin, std::endl, std::string;

string FILE_CHOICE = "../dataset/relA_mock.txt";

int main() {

    cout << "****************************************" << endl;
    cout << "********* DNA Sequence Finder **********" << endl;
    cout << "****************************************" << endl;

    char menuChoice = 'Y';
    do {
        // open the file
        std::ifstream file(FILE_CHOICE);
        if (!file.is_open()) {
            cout << "Error: file not found" << endl;
            return 1;
        }

        //get the sequence to find
        string sequenceToFind = getSequenceToFind();

        std::atomic<int> matchCount(0);
        std::vector<std::future<std::string>> futures;

        string line;
        while (std::getline(file, line)) {
            futures.push_back(std::async(std::launch::async, [line, &sequenceToFind, &matchCount]() mutable {
                size_t pos = 0;
                int matchCountPerLine = 0;
                std::vector<size_t> matchPositions;
                while ((pos = line.find(sequenceToFind, pos)) != string::npos) {
                    matchCountPerLine++;
                    matchPositions.push_back(pos);
                    pos += sequenceToFind.length();
                }
                matchCount += matchCountPerLine;
                for (auto it = matchPositions.rbegin(); it != matchPositions.rend(); ++it) {
                    line.insert(*it, "\033[33m");
                    line.insert(*it + sequenceToFind.length() + 5, "\033[0m");
                }
                return line;
            }));
        }

        std::vector<std::string> lines;
        for (auto& future : futures) {
            lines.push_back(future.get());
        }
        for (const auto& line : lines) {
            cout << line << endl;
        }

        cout << "Total sequences found: " << matchCount << endl;
        //go again?
        menuChoice = getGoAgainChoice();
    } while (menuChoice != 'N' && menuChoice != 'n');
    cout << "Goodbye!" << endl;

    return 0;
}