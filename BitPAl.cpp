#include "BitPAL.h"
#include <iostream> 
#include <cstdlib>
#include <fstream> // used to read files
#include <vector>
#include <cstdint>

using namespace std;

// Hard coded match/mismatch score and gap penalty, by consequence the minimum and maximum value for delta_H/V are respectively -3 and 4
int M = 1;
int I = -1;
int G = -3; 


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

// deltaVmax_shift = deltaVpos7_shift (algorithm 1 seen during the lecture)
uint64_t create_deltaVmax_shift(uint64_t deltaHmin, uint64_t match){
    uint64_t Init_pos_max = 0;
    uint64_t deltaVmax_shift = 0;
    Init_pos_max = deltaHmin & match;
    deltaVmax_shift = ((Init_pos_max + deltaHmin)^deltaHmin)^Init_pos_max;
    return deltaVmax_shift;
}

// deltaVhigh(algorithm 2 seen during the lecture)
uint64_t create_deltaVhigh(uint64_t bit_vectors_delta_H){} //Afwerken

// deltaVlow (algorithm 3 seen during the lecture)
uint64_t create_deltaVlow(){} //Afwerken

// deltaVmin (algorithm 4 seen during the lecture)
uint64_t create_deltaVmin(vector<uint64_t> bit_vectors_delta_V){
    uint64_t deltaHmin_shift = 0;
    uint64_t all_ones = ~0ULL;
    uint64_t or_result = 0;
    for (int i = G + 1; i <= M-G; i++) {
        or_result |= bit_vectors_delta_V[i];
    }
    deltaHmin_shift = all_ones^(or_result);
}

int main(){
    vector<string> sequences;
    read_fasta_file("global_2.fasta", sequences);

    vector<uint64_t> match_vectors(4,0); 
    create_match_vectors(match_vectors, sequences);

    vector<uint64_t> bit_vectors_delta_H(M-G-G,0);
    vector<uint64_t> bit_vectors_delta_V(M-G-G,0);
}