#include "DNAFinder.h"
#include <iostream>
#include <fstream>
#include <string>

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

        string line;
        int matchCount = 0;
        while (std::getline(file, line)) {
            size_t pos = 0;
            int matchCountPerLine = 0;
            while ((pos = line.find(sequenceToFind, pos)) != string::npos) {
                matchCountPerLine++;
                line.insert(pos, "\033[33m");
                pos += sequenceToFind.length() + 5;
                line.insert(pos, "\033[0m");
                pos += 4;
            }
            matchCount += matchCountPerLine;
            cout << line << endl;
        }
        cout << "Total sequences found: " << matchCount << endl;
        //go again?
        menuChoice = getGoAgainChoice();
    } while (menuChoice != 'N' && menuChoice != 'n');
    cout << "Goodbye!" << endl;
   

    return 0;
}