1. Compile the program using g++: g++ -O2 BitPAl.cpp -o bitpal
2. Run the executable: ./bitpal input.fasta
3. The output (the sequences to be aligned and the global alignment score) will appear on the screen. 

Remark: the code is not error free. There is a bug in the calculation of deltaVlow resulting in a bad alignment score for some sequences. (e.g. sequences resulting in a good alignment score: sequence 1: ACCGA, sequence 2:ACCGC) 