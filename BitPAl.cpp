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

int mid = I-G;

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

// deltaVhigh(algorithm 2 seen during the lecture) (mid+1 to max-1)
uint64_t create_deltaVhigh(vector<uint64_t> bit_vectors_delta_H,vector<uint64_t>& bit_vectors_delta_V_shift, uint64_t match){
    uint64_t remaindeltaHmin = 0;
    remaindeltaHmin = bit_vectors_delta_H[0]^(bit_vectors_delta_V_shift.back() >> 1); //klopt dit?

    uint64_t deltaVmax_shift_or_match = bit_vectors_delta_V_shift.back() | match;

    for (int i = bit_vectors_delta_H.size() - 2; i > mid - G; i--){ // "V=max - 1 to "V=mid+1 (correct)
        uint64_t initpos = bit_vectors_delta_H[M-G-i] & deltaVmax_shift_or_match; // 1st term = max - the value of interest
        for (int j = bit_vectors_delta_H.size() - 3; j > i; i--){
            uint64_t deltaVpos_shift_not_match = bit_vectors_delta_V_shift[j] & ~match; //1st term must be the shifted one!!
            initpos = initpos | (bit_vectors_delta_H[M-G-j] & deltaVpos_shift_not_match);
        }
        uint64_t deltaVshift = ((initpos << 1)+remaindeltaHmin)^remaindeltaHmin;
        uint64_t deltaVshift_not_match = deltaVshift & ~match;
        bit_vectors_delta_V_shift[i] =  deltaVshift_not_match;
    }

} //Controleren

// deltaVlow (algorithm 3 seen during the lecture)
uint64_t create_deltaVlow(vector<uint64_t> bit_vectors_delta_V_shift, uint64_t match){
    uint64_t deltaVmax_shift_or_match = bit_vectors_delta_V_shift.back() | match;
    uint64_t or_result = 0;
    for (int i = mid - G - 1; i > 0; i--){ // "V = mid to "V = min + 1
        for (int i = bit_vectors_delta_V_shift.size() - 2; i > mid - G - 1; i--){
            or_result = or_result | bit_vectors_delta_V_shift[i];
        }
        uint64_t deltaVnotmaxtomidplusoneshiftormatch = ~(deltaVmax_shift_or_match|or_result);
    }
}// afwerken

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

    vector<uint64_t> bit_vectors_delta_H(M-G-G+1,0);
    vector<uint64_t> bit_vectors_delta_V(M-G-G+1,0);
}