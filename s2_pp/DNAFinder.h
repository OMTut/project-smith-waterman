
#ifndef DNA_FINDER_H
#define DNA_FINDER_H

#include <iostream>
#include <fstream>
#include <string>

using std::cout, std::cin, std::endl, std::string;

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