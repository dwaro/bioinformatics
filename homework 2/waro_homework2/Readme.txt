The source code to align sequences lives in align_proteins.py .

The program can be run using the following command:
python3 align_proteins.py seq1_file.fna seq2_file.fna [matches.txt]

e.g. python3 align_proteins.py SARS-CoV-1_N_protein.fna SARS-CoV-2_N_protein.fna

The permute.py file is the code for doing the permutations, but is not needed to
do the alignments. Seq1_file.fna and seq2_file.fna are required to run the program,
while the matches.txt file optional.

The output of the alignments is written to output.txt **Note this file does overwrite
itself with every execution of the program.
