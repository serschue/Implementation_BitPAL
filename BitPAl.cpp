#include "BitPAL.h"
#include <iostream> 
#include <cstdlib>
#include <fstream> // used to read files
#include <vector>

using namespace std;

void read_fasta_file(const string& file_name, vector<string>& sequences){
    ifstream file(file_name.c_str());
    if (!file){
        throw runtime_error("The file: " + file_name + "could not be found");
    }

    string line;
    while (getline(file,line)) {
        if (line.empty()){
            continue;
        }
        if (line.front() == '>'){
            sequences.push_back(string());
            continue;
        }

        sequences.back().append(line);
        cout << "Sequence: " << sequences.back() << "\n";
    }

    file.close();
}

int main(){
    vector<string> sequences;
    read_fasta_file("global.fasta", sequences);
}