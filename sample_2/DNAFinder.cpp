#include "DNAFinder.h"
#include <iostream>
#include <fstream>
#include <string>

using std::cout, std::cin, std::endl, std::string;

string FILE_CHOICE;

int main() {

    cout << "****************************************" << endl;
    cout << "********* DNA Sequence Finder **********" << endl;
    cout << "****************************************" << endl;

    FILE_CHOICE = chooseFileSize();

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

        //find sequence
        bool found = false;
        int matchCount = 0;
        int lineCount = 1;
        string line;
        while (getline(file, line)) {
            size_t pos = 0;
            string highlightedLine = line;
            int matchCountPerLine = 0;
            while ((pos = highlightedLine.find(sequenceToFind, pos)) != string::npos) {
                matchCountPerLine++;
                highlightedLine.insert(pos, "\033[33m");
                pos += 5;
                highlightedLine.insert(pos + sequenceToFind.length(), "\033[0m");
                pos += sequenceToFind.length() + 5;
            }

            // print the line if a match is found
            if (matchCountPerLine > 0) {
                matchCount++;
                cout << "Sequence found on line " << lineCount << ": " << highlightedLine << endl;
                found = true;
            }
            lineCount++;
        }

        // print out the total sequences with a match
        totalSequencesFound(found, matchCount);


        //go again?
        menuChoice = getGoAgainChoice();
    } while (menuChoice != 'N' && menuChoice != 'n');
    cout << "Goodbye!" << endl;
   

    return 0;
}