#include "BitPAL.h"
#include <iostream> 
#include <cstdlib>
#include <fstream> // used to read files
#include <vector>
#include <cstdint>

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

void create_match_vectors(vector<uint64_t>& match_vectors, vector<string> sequences){
    string horizontal_sequence = sequences[0];
    for (int i = sequences[0].length() - 1; i >= 0; i--){
        if (horizontal_sequence[i] == 'A'){ // match vector A
            match_vectors[0] |= (1ULL << i);
        }
        if (horizontal_sequence[i] == 'C'){ // match vector C
            match_vectors[1] |= (1ULL << i);
        }
        if (horizontal_sequence[i] == 'G'){ // match vector G
            match_vectors[2] |= (1ULL << i);
        }
        if (horizontal_sequence[i] == 'T'){ // match vector T
            match_vectors[3] |= (1ULL << i);
        }
    }
    cout << match_vectors[2];
}


int main(){
    vector<string> sequences;
    read_fasta_file("global_2.fasta", sequences);

    vector<uint64_t> match_vectors(4,0); 
    create_match_vectors(match_vectors, sequences);
}