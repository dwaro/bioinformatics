This program is written in python 3.

Executing the program from the command line:
>>> python3.8 main.py path_to_hw3_file.fna

Outputs:
- edges.txt
- tree.txt
- genetic-distances.txt
- bootstrap.txt

To build the tree images for steps 4 and 5, the following R scripts were used provided by Professor Knight:
- (step 4): Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt
- (step 5): Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt bootstrap.txt