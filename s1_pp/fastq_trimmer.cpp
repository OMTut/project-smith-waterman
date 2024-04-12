#include <iostream>
#include <fstream>
#include <string>

int main() {

   using std::cout, std::cin, std::endl, std::string;

   std::ifstream file("../dataset/ERR000001_FASTQ.txt");
   std::ofstream output("../dataset/ERR000001_FASTQ_FULL.txt");

   cout << "Running..." << endl;

   if (!file.is_open()) {
      cout << "Error: file not found" << endl;
      return 1;
   }

   string line;
   bool printNextLine = false;
   int count = 0;
   while (getline(file, line)) {
      if (printNextLine) {
         output << line << endl;
         printNextLine = false;
         count++;
      } else if (!line.empty() && line[0] == '@') {
         printNextLine = true;
      }
   }
   output.close();
   return 0;
}