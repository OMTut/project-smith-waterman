
#ifndef DNA_FINDER_H
#define DNA_FINDER_H

#include <iostream>
#include <fstream>
#include <string>

using std::cout, std::cin, std::endl, std::string;


/*****************************
 * chooseFileSize
 * return string
 * Desc: returns the file path as a string
 *       1 for Test File (100 sequences) or 2 for Full File
******************************/
string chooseFileSize() {
    int fileChoice;
    do {
        cout << "Enter 1 for Test File (100 sequences) or 2 for Full File (100000 sequences): ";
        cin >> fileChoice;
        if (fileChoice != 1 && fileChoice != 2) {
            cout << "Invalid choice. Please try again." << endl;
        }
    } while (fileChoice != 1 && fileChoice != 2);
    if (fileChoice == 1) {
        return "../dataset/ERR000001_FASTQ_100.txt";
    } else {
        return "../dataset/ERR000001_FASTQ_FULL.txt";
    }
}


/*****************************
 * getSequenceToFind
 * return string
 * Desc: Prompts the user to enter a DNA sequence to find
******************************/
string getSequenceToFind() {
   string sequenceToFind;
   cout << "Enter a DNA sequence to find: ";
   cin >> sequenceToFind;
   return sequenceToFind;
}

void totalSequencesFound(bool found, int matchCount) {
    if (!found) {
        cout << "Sequence not found" << endl;
    } else {
        cout << "Total sequences found: " << matchCount << endl;
    }
}

/*****************************
 * getGoAgainChoice
 * return char
 * Desc: Prompts the user to enter Y or N to search for another sequence
******************************/
char getGoAgainChoice() {
   char goAgainChoice;
   cout << "Would you like to search for another sequence? (Y/N): ";
   cin >> goAgainChoice;
   return goAgainChoice;
}



#endif // DNA_FINDER_H