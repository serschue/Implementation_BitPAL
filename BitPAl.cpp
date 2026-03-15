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
int min = G;

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

void printBitPattern(uint64_t value) {
    for (int i = 63; i >= 0; i--) {
        std::cout << ((value >> i) & 1);
    }
    std::cout << std::endl;
}

void create_match_vectors(vector<uint64_t>& match_vectors, vector<string> sequences){
    string horizontal_sequence = sequences[0];
    for (int i = sequences[0].length() - 1; i >= 0; i--){
        uint64_t mask = (1ULL << i);
        if (horizontal_sequence[i] == 'A'){ // match vector A
            match_vectors[0] |= mask;
        }
        if (horizontal_sequence[i] == 'C'){ // match vector C
            match_vectors[1] |= mask;
        }
        if (horizontal_sequence[i] == 'G'){ // match vector G
            match_vectors[2] |= mask;
        }
        if (horizontal_sequence[i] == 'T'){ // match vector T
            match_vectors[3] |= mask;
        }
        // debug
        //printBitPattern(match_vectors[0]);
        //printBitPattern(match_vectors[1]);
        //printBitPattern(match_vectors[2]);
        //printBitPattern(match_vectors[3]);    
    }
}

// deltaVmax_shift (algorithm 1 seen during the lecture)
void create_deltaVmax_shift(vector<uint64_t>& bit_vectors_delta_V_shift_curr,uint64_t deltaHmin, uint64_t match){
    uint64_t Init_pos_max = 0;
    uint64_t deltaVmax_shift = 0;
    Init_pos_max = deltaHmin & match;
    deltaVmax_shift = ((Init_pos_max + deltaHmin)^deltaHmin)^Init_pos_max;
    bit_vectors_delta_V_shift_curr.back() = deltaVmax_shift;
}

// deltaHmax_shift
void create_deltaHmax(vector<uint64_t>& bit_vectors_delta_H_curr,vector<uint64_t>& bit_vectors_delta_H_prev, vector<uint64_t>& bit_vectors_delta_V_shift_curr,const vector<string>& sequences){
    uint64_t deltaHmax = 0;
    deltaHmax = (bit_vectors_delta_V_shift_curr.back() & bit_vectors_delta_H_prev[0]);
    uint64_t all_ones = (1ULL << sequences[0].size()+1) - 1;
    printBitPattern(all_ones);
    printBitPattern(bit_vectors_delta_V_shift_curr.back());
    printBitPattern(bit_vectors_delta_H_prev[0]);
    bit_vectors_delta_H_curr.back() = all_ones^deltaHmax;
}

// deltaVhigh (algorithm 2 seen during the lecture) (mid+1 to max-1)
void create_deltaVhigh(vector<uint64_t> bit_vectors_delta_H,vector<uint64_t>& bit_vectors_delta_V_shift, uint64_t match){
    uint64_t remaindeltaHmin = bit_vectors_delta_H[0]^(bit_vectors_delta_V_shift.back() >> 1);

    uint64_t deltaVmax_shift_or_match = bit_vectors_delta_V_shift.back() | match;

    for (int i = bit_vectors_delta_H.size() - 2; i > mid - G; i--){ // "V = max - 1 to "V=mid+1 (correct)
        uint64_t initpos = bit_vectors_delta_H[M-G-i] & deltaVmax_shift_or_match; // 1st term = max - the value of interest (correct)
        // Is not executed due to the choice of M/I and G
        for (int j = bit_vectors_delta_H.size() - 2; j > i; j--){ // max-1 to the high value you are calculating
            uint64_t deltaVpos_shift_not_match = bit_vectors_delta_V_shift[j] & ~match; //1st term must be the shifted one!!
            initpos = initpos | (bit_vectors_delta_H[M-G-j] & deltaVpos_shift_not_match);
        }
        uint64_t deltaVshift = ((initpos << 1)+remaindeltaHmin)^remaindeltaHmin;
        uint64_t deltaVshift_not_match = deltaVshift & ~match;
        bit_vectors_delta_V_shift[i] =  deltaVshift_not_match;
    }
}

// deltaHhigh (mid+1 to max-1)
void create_deltaHhigh_probeer(vector<uint64_t> bit_vectors_delta_H,vector<uint64_t>& bit_vectors_delta_V_shift_curr, uint64_t match){ 
    int k = 0;
    for (int i = M-G-1; i > mid; i--){ // going from max - 1 = 3 to mid
        k = 1;
        uint64_t or_result = match & (bit_vectors_delta_H[k] << 1);
        for (int j = M - G - 1; j > i - 1; j--){ 
            or_result = or_result | ((bit_vectors_delta_H[k] << 1) &(bit_vectors_delta_V_shift_curr[j-G] & ~match));
            k = k - 1;
        }
        bit_vectors_delta_H[i-G] =  or_result;
    }
}

// deltaHhigh (mid+1 to max-1)
void create_deltaHhigh(vector<uint64_t>& bit_vectors_delta_H_prev,vector<uint64_t>& bit_vectors_delta_H_curr, vector<uint64_t>& bit_vectors_delta_V_shift_curr, uint64_t match){ 
    int k = 0;
    for (int i = M-G-1; i > mid; i--){ // going from max - 1 = 3 to mid
        k = 1;
        uint64_t or_result = 0;
        for (int j = M - G - 1; j > i - 1; j--){ 
            or_result = or_result | ((bit_vectors_delta_H_prev[k]) & (bit_vectors_delta_V_shift_curr[j-G]));
            k = k - 1;
        }
        bit_vectors_delta_H_curr[i-G] =  or_result;
    }
}

// deltaVlow (algorithm 3 seen during the lecture)
void create_deltaVlow_eerdere_versie(vector<uint64_t>& bit_vectors_delta_V_shift,const vector<uint64_t>& bit_vectors_delta_H, uint64_t match){
    uint64_t deltaVmax_shift_or_match = bit_vectors_delta_V_shift.back() | match;
    uint64_t or_result = deltaVmax_shift_or_match;
    for (int i = mid - G; i > 0; i--){ // "V = mid to "V = min + 1
        uint64_t deltaVlow_shift = 0;
        uint64_t a = 0;
        for (int j = bit_vectors_delta_V_shift.size() - 2; j > i; j--){ // 1 + the value of the bitvector you are calculating
            or_result = or_result | bit_vectors_delta_V_shift[j];
            deltaVlow_shift = deltaVlow_shift | (bit_vectors_delta_H[j-a+2*G] & (bit_vectors_delta_V_shift[j]&~match)); //Nog een fout in bit_vectors_deltaH
        }
        uint64_t deltaVnotmaxtomidplusoneshiftormatch = ~(or_result);
        bit_vectors_delta_V_shift[i] = (deltaVlow_shift|(bit_vectors_delta_H[0] & deltaVnotmaxtomidplusoneshiftormatch)) << 1;

    }
}

// deltaVlow (algorithm 3 seen during the lecture)
void create_deltaVlow(vector<uint64_t>& bit_vectors_delta_V_shift,const vector<uint64_t>& bit_vectors_delta_H, uint64_t match){
    uint64_t deltaVmax_shift_or_match = bit_vectors_delta_V_shift.back() | match;
    uint64_t or_result = deltaVmax_shift_or_match;
    for (int i = mid - G; i > 0; i--){ // "V = mid to "V = min + 1 (correct)
        uint64_t deltaVlow_shift = 0;
        uint64_t a = 0;
        for (int j = bit_vectors_delta_V_shift.size() - 2; j > i; j--){ // 1 + the value of the bitvector you are calculating
            or_result = or_result | bit_vectors_delta_V_shift[j];
        }
        uint64_t k = 0;
        int t = i + k;
        while(t < M-G){
            deltaVlow_shift = deltaVlow_shift | (bit_vectors_delta_H[k] & (bit_vectors_delta_V_shift[i+k-G]&~match));
            k = k+1;
            t = i + k;
        }
        uint64_t deltaVnotmaxtomidplusoneshiftormatch = ~(or_result);
        bit_vectors_delta_V_shift[i] = (deltaVlow_shift|(bit_vectors_delta_H[0] & deltaVnotmaxtomidplusoneshiftormatch)) << 1;
    }
}

// deltaHlow
void create_deltaHlow_probeer(vector<uint64_t>& bit_vectors_delta_V, const vector<uint64_t>& bit_vectors_delta_H, uint64_t match){
    for (int i = mid; i > 0; i--){ // going from mid = 2 to min + 1
        uint64_t or_result = match & (bit_vectors_delta_H[M-G - i] << 1); // max - the value you want to compute
        int k = 1;
        for (int j = 0; j < bit_vectors_delta_H.size(); j++){
            if (j < mid + 1){
                or_result = or_result | ((bit_vectors_delta_H[0] << 1) & (bit_vectors_delta_V[j] & ~match));
            }
            else{
                or_result = or_result | ((bit_vectors_delta_H[k] << 1) & (bit_vectors_delta_V[j] & ~match));
                k = k + 1;
            }
        }
    }
}

// deltaHlow
void create_deltaHlow(vector<uint64_t>& bit_vectors_delta_V_shift_curr, const vector<uint64_t>& bit_vectors_delta_H_prev, vector<uint64_t>& bit_vectors_delta_H_curr, uint64_t match){
    for (int i = mid; i > 0; i--){ // going from mid = 2 to min + 1
        uint64_t or_result = 0;
        int k = 1;
        for (int j = 0; j < bit_vectors_delta_V_shift_curr.size(); j++){
            if (j < mid + 1 - G){
                or_result = or_result | ((bit_vectors_delta_H_prev[0]) & (bit_vectors_delta_V_shift_curr[j]));
            }
            else{
                or_result = or_result | ((bit_vectors_delta_H_prev[k]) & (bit_vectors_delta_V_shift_curr[j]));
                k = k + 1;
            }
        }
    }
}

// deltaVmin (algorithm 4 seen during the lecture)
void create_deltaVmin(const vector<string>& sequences, vector<uint64_t>& bit_vectors_delta_V_shift){
    uint64_t deltaVmin_shift = 0;
    uint64_t all_ones = (1ULL << sequences[0].size()+1) - 1;
    uint64_t or_result = 0;
    for (int i = 1; i < bit_vectors_delta_V_shift.size(); i++) {
        or_result |= bit_vectors_delta_V_shift[i];
    }
    bit_vectors_delta_V_shift[0] = all_ones^(or_result);
}

// deltaHmin
void create_deltaHmin(const vector<string>& sequences, vector<uint64_t>& bit_vectors_delta_H_curr){
    uint64_t deltaHmin = 0;
    uint64_t all_ones = (1ULL << sequences[0].size()) - 1;
    uint64_t or_result = 0;
    for (int i = 1; i < bit_vectors_delta_H_curr.size(); i++) {
        or_result |= bit_vectors_delta_H_curr[i];
    }
    printBitPattern(bit_vectors_delta_H_curr[7]);
    bit_vectors_delta_H_curr[0] = all_ones^(or_result);
}

int calculate_global_alignment_score(const vector<string>& sequences, const vector<uint64_t>& bit_vectors_delta_H){
    int global_alignment_score = sequences[1].size()*G;
    for (int i = 0; i < sequences[0].size(); i++){
        for (int j = 0; j < bit_vectors_delta_H.size(); j++){
            if (((bit_vectors_delta_H[j] >> i) & 1) == 1){
                global_alignment_score += j;
            }
        }
    }
    return global_alignment_score;
}

void printDeltaValues(const vector<uint64_t>& bit_vectors) {
    for (int i = 0; i < 64; i++) {
        uint64_t mask = 1ULL << i;

        for (int j = 0; j < 8; j++) {
            if (bit_vectors[j] & mask) {
                int value = j - 3;
                printf("%2d ", value);
                break;
            }
        }
    }
    printf("\n");
}

int main(){
    vector<string> sequences;
    read_fasta_file("global_2.fasta", sequences);

    vector<uint64_t> match_vectors(4,0); 
    create_match_vectors(match_vectors, sequences); 

    vector<uint64_t> bit_vectors_delta_H_prev(M-G-G+1,0);
    vector<uint64_t> bit_vectors_delta_H_curr(M-G-G+1,0);
    vector<uint64_t> bit_vectors_delta_V_shift_prev(M-G-G+1,0);
    vector<uint64_t> bit_vectors_delta_V_shift_curr(M-G-G+1,0);

    string horizontal_sequence = sequences[0];
    string vertical_sequence = sequences[1];

    // Initalisation: set the bitvectors for delta_H_previous
    bit_vectors_delta_H_prev[0] = (1ULL << horizontal_sequence.size()) - 1;

    // Iteration over the vertical sequence
    for (int i = 3; i < 5; i++){//vertical_sequence.size(); i++){

        // Set bit_vectors of delta_H_curr to zero
        std::fill(bit_vectors_delta_H_curr.begin(), bit_vectors_delta_H_curr.end(), 0);
        std::fill(bit_vectors_delta_V_shift_curr.begin(), bit_vectors_delta_V_shift_curr.end(), 0);

        // Set the correct match vector
        uint64_t match = 0;
        if (vertical_sequence[i]=='A'){
            match = match_vectors[0];
        }
        else if (vertical_sequence[i]=='C'){
            match = match_vectors[1];
        }
        else if (vertical_sequence[i]=='G'){
            match = match_vectors[2];
        }
        else{
            match = match_vectors[3];
        }
        cout << i << endl;

        create_deltaVmax_shift(bit_vectors_delta_V_shift_curr,bit_vectors_delta_H_prev[0], match);
        // Until here correct
        create_deltaVhigh(bit_vectors_delta_H_prev, bit_vectors_delta_V_shift_curr,match);
        create_deltaVlow(bit_vectors_delta_V_shift_curr,bit_vectors_delta_H_prev,match);
        cout << "DV4" << endl;
        for (int i = 0; i < 8; i++){
            printBitPattern(bit_vectors_delta_V_shift_curr[i]);
        }
        cout << "--------" << endl;
        printDeltaValues(bit_vectors_delta_V_shift_curr);
        create_deltaVmin(sequences,bit_vectors_delta_V_shift_curr);
        // I think until here correct
        printDeltaValues(bit_vectors_delta_V_shift_curr);

        create_deltaHmax(bit_vectors_delta_H_curr,bit_vectors_delta_H_prev,bit_vectors_delta_V_shift_curr,sequences);
        create_deltaHhigh(bit_vectors_delta_H_prev,bit_vectors_delta_H_curr, bit_vectors_delta_V_shift_curr,match);          
        create_deltaHlow(bit_vectors_delta_V_shift_curr,bit_vectors_delta_V_shift_prev,bit_vectors_delta_H_curr,match);
        create_deltaHmin(sequences,bit_vectors_delta_H_curr); 
        
        printDeltaValues(bit_vectors_delta_H_curr);
        bit_vectors_delta_H_prev = bit_vectors_delta_H_curr;
        bit_vectors_delta_V_shift_prev = bit_vectors_delta_V_shift_curr;
    }
    //int global_alignment_score = calculate_global_alignment_score(sequences, bit_vectors_delta_H);
}