/******************************************************************
 * <mutex>
 * a library for mutexes which are used to protect shared data from
 * being accessed by multiple threads at the same time.
 * 
 * <future>
 * a library for asyncronous programming. It allows you to access the
 * results of asyncronous functions that might not have completed yet.
 * When the async launches, it returns a future object that can check if
 * a result is ready. It will block until the result is ready.
*/

#include "DNAFinder.h"
#include <iostream>
#include <fstream>
#include <string>
#include <future>
#include <vector>
#include <atomic>
#include <algorithm>
#include <thread>
#include <mutex>

using std::cout, std::cin, std::endl, std::string;

string FILE_CHOICE;

int main() {

    cout << "****************************************" << endl;
    cout << "********* DNA Sequence Finder **********" << endl;
    cout << "****************************************" << endl;

    

    char menuChoice = 'Y';
    do {
        // open the file
        FILE_CHOICE = chooseFileSize();
        std::ifstream file(FILE_CHOICE);
        if (!file.is_open()) {
            cout << "Error: file not found" << endl;
            return 1;
        }

        //get the sequence to find
        string sequenceToFind = getSequenceToFind();

        std::atomic<int> matchCount(0);
        std::atomic<int> lineCount(0);

        //chunksize: the number of lines each thread will process
        const int chunkSize = 250000;
        std::vector<std::string> lines(chunkSize);
        std::vector<std::thread> threads;

        //allLines stores pairs. the line number found and the line itself
        std::vector<std::pair<int, std::string>> allLines;
        
        //allLinesMutex is used to protect the allLines vector from
        //being accessed by multiple threads at the same time
        std::mutex allLinesMutex;
        // holds a tuple of the match count, line count, and the lines that contain the sequence
        // is what will build allLines
        std::vector<std::future<std::tuple<int, int, std::vector<std::pair<int, std::string>>>>> futures;

        while (true) {
            int i;
            for (i = 0; i < chunkSize && std::getline(file, lines[i]); ++i);
            if (i == 0) break;

            /**********************************
             * async launches and creates a new thread
            */
            futures.push_back(std::async(std::launch::async, [lines, &sequenceToFind, i]() mutable {
                int localMatchCount = 0;
                int localLineCount = 0;
                std::vector<std::pair<int, std::string>> result;
                for (int j = 0; j < i; ++j) {
                    localLineCount++;
                    std::string& line = lines[j];
                    size_t pos = 0;
                    int matchCountPerLine = 0;
                    while ((pos = line.find(sequenceToFind, pos)) != std::string::npos) {
                        matchCountPerLine++;
                        line.insert(pos, "\033[33m");
                        pos += 5;
                        line.insert(pos + sequenceToFind.length(), "\033[0m");
                        pos += sequenceToFind.length() + 5;
                    }
                    if (matchCountPerLine > 0) {
                        localMatchCount++;
                        result.push_back(std::make_pair(localLineCount, line));
                    }
                }
                return std::make_tuple(localMatchCount, localLineCount, result);
            }));
        }

        int totalMatchCount = 0;
        int totalLineCount = 0;
        for (auto& future : futures) {
            int localMatchCount, localLineCount;
            std::vector<std::pair<int, std::string>> result;
            std::tie(localMatchCount, localLineCount, result) = future.get();
            totalMatchCount += localMatchCount;
            totalLineCount += localLineCount;
            std::lock_guard<std::mutex> guard(allLinesMutex);
            allLines.insert(allLines.end(), result.begin(), result.end());
        }

        std::sort(allLines.begin(), allLines.end());

        // for (const auto& result : allLines) {
        //     cout << "Sequence found on line " << result.first << ": " << result.second << endl;
        // }
        cout << "Total sequences found: " << totalMatchCount << endl;

        //go again?
        menuChoice = getGoAgainChoice();
    } while (menuChoice != 'N' && menuChoice != 'n');
    cout << "Goodbye!" << endl;

    return 0;
}