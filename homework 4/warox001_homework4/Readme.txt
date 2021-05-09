The code file (main.py) calculates/plots conservation rates and calculates/plots variable regions
for a given aligned sequence file with primers. Using a set of parameters (minimum region size, 
maximum region size, and number of variable regions), the program calculates the variances of all
possible sections of the gene positions and chooses the variable regions with the most variance.

The program can be run via Python3 with the following command:
>>> Python3 main.py path_to_seqs_with_primers.fna

The outputs are:
- solution-problem-1.txt (file of conservation rates for each gene position)
- solution-problem-2.pdf (plot of smoothed conservation rates for each gene position)
- solution-problem-3.txt (file of start and end positions for each variable region)
- solution-problem-4.pdf (plot of variable regions on top of conservation rates)