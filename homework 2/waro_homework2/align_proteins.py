import sys, re

def main(seq1_file, seq2_file, matches_file=None):
    """
    Description:
        This function is responsible for orchestrator the flow of the program.

    Parameters:
        - seq1_file (str): file path of first sequence file
        - seq2_file (str): file path of second sequence file
        - matches_file (str)(optional): the file containing where sequences 1 &
            2 are known to align. Default to None.
    """
    seq1 = open_file(seq1_file)
    seq2 = open_file(seq2_file)

    if matches_file is None:
        # no matches file included with input. Just run needleman-wunsch on all
        # of sequences 1 and 2
        score, aligned1, aligned2 = needleman_wunsch(seq1, seq2)
    else:
        # matches file included, only run needleman-wunsch on non-anchored
        # sections of the sequences
        score, aligned1, aligned2 = anchored(seq1, seq2, matches_file)

    print(f"Score: {score}")
    save_file(score, aligned1, aligned2, seq1_file, seq2_file) # save file to output.txt
    print(f"Aligned sequences saved to output.txt")

def needleman_wunsch(seq1, seq2):
    """
    Description:
        This function is responsible for performing the needleman-wunsch
        alignment algorithm on the two input sequences.

    Parameters:
        - seq1 (str): sequence 1 of 2 to be aligned
        - seq2 (str): sequence 2 of 2 to be aligned

    Returns tuple(score, aligned1, aligned2):
        - score (int): the alignment score of aligning seq1 and seq2
        - aligned1 (str): seq1 adjusted to be aligned with seq2
        - aligned2 (str): seq2 adjusted to be aligned with seq1
    """
    width = len(seq1) # width of the matrix
    height = len(seq2) # height of the matrix
    matrix = initialize_matrix(width, height) # matrix containing alignment scores
    gap = -2 # gap penalty
    # ignore the affine/linar gap penalty. i.e. all gaps are worth -2

    # populate the matrix
    for row in range(1, height + 1):
        for col in range(1, width + 1):
            # diagonally, will it be a match(+1) or a mismatch(-3)
            match = 1 if seq1[col - 1] == seq2[row - 1] else -3

            # assign the max of diagonal, or gaps from above/left to next cell
            matrix[row][col] = max(
                matrix[row-1][col] + gap,
                matrix[row][col-1] + gap,
                matrix[row-1][col-1] + match
            )

    score = matrix[height][width] # grab the score from the last cell

    # move back through the matrix, aligning the sequences
    aligned1, aligned2 = align_sequences(matrix, seq1, seq2)

    return (score, aligned1, aligned2)

def align_sequences(matrix, seq1, seq2):
    """
    Description:
        This function steps backwards through the alignment matrix, choosing the
        optimal path from the last cell (bottom-right) to the starting cell.
        This results with the optimal alignment between seq1 and seq2.

    Parameters:
        - matrix list(list): 2D alignment matrix populated with scores for
            aligning seq1 and seq2
        - seq1 (str): sequence 1 of 2 being aligned
        - seq2 (str): sequence 2 of 2 being aligned

    Returns tuple(aligned1, aligned2):
        - aligned1 (str): seq1 adjusted to be aligned with seq2
        - aligned2 (str): seq2 adjusted to be aligned with seq1
    """
    aligned1 = ""
    aligned2 = ""

    # start from last row and column
    row = len(matrix) - 1
    col = len(matrix[0]) - 1

    # loop from the final position until we're back to the initial position.
    # we take a greedy path back, seeing if it looks like we came from above or
    # from the left, and if not, then we came from the diagonal
    while row > 0 or col > 0:
        current = matrix[row][col] # value of our current spot in the matrix

        # we came from above, add gap to seq1
        if row > 0 and current == matrix[row-1][col] - 2:
            aligned1 = '-' + aligned1
            aligned2 = seq2[row-1] + aligned2
            row -= 1

        # we came from the left, add gap to seq2
        elif col > 0 and current == matrix[row][col-1] - 2:
            aligned1 = seq1[col-1] + aligned1
            aligned2 = '-' + aligned2
            col -= 1

        # we came from diagonal
        else:
            aligned1 = seq1[col-1] + aligned1
            aligned2 = seq2[row-1] + aligned2
            row -= 1
            col -= 1

    return (aligned1, aligned2)

def anchored(seq1, seq2, matches_file):
    """
    Description:
        This function is responsible for performing the anchored alignment of
        seq1 and seq2.

    Parameters:
        - seq1 (str): sequence 1 of 2 to be aligned
        - seq2 (str): sequence 2 of 2 to be aligned
        - matches_file (str): the file path for the matches file containing
            where seq1 and seq2 are known to align.

    Returns tuple(score, aligned1, aligned2):
        - score (int): the alignment score of aligning seq1 and seq2
        - aligned1 (str): seq1 adjusted to be aligned with seq2
        - aligned2 (str): seq2 adjusted to be aligned with seq1
    """
    # extract anchored locations for the different CoV proteins
    cov1, cov2 = open_matches_file(matches_file)

    aligned1, aligned2 = "", ""

    # align up to the first anchored alignments
    score, aligned1, aligned2 = needleman_wunsch(seq1[0:cov1[0][0]-1], seq2[0:cov2[0][0]-1])

    # iterate over the anchored alignments, first simply adding the anchored
    # alignment sections, and then performing needleman-wunsch for the sections
    # between achored sections.
    for i in range(len(cov1)):
        # add anchored alignment areas
        aligned1 += seq1[cov1[i][0]-1:cov1[i][1]]
        aligned2 += seq2[cov2[i][0]-1:cov2[i][1]]

        if i < len(cov1) - 1:
            # if there are more anchored sections, perform needleman-wunsch
            # between this anchored section and the next one.
            score, a1, a2 = needleman_wunsch(seq1[cov1[i][1]:cov1[i+1][0]-1], seq2[cov2[i][1]:cov2[i+1][0]-1])
        else:
            # no more anchored sections. Do needleman-wunsch on the remaining
            # parts of seq1 and seq2
            score, a1, a2 = needleman_wunsch(seq1[cov1[i][1]:], seq2[cov2[i][1]:])

        # add non-anchored sections
        aligned1 += a1
        aligned2 += a2

    # calculate the score of the alignment
    score = 0
    for i in range(len(aligned1)):
        aa1 = aligned1[i] # seq1 char at position i
        aa2 = aligned2[i] # seq2 cahr at position i

        if aa1 == '-' or aa2 == '-':
            score -= 2 # gap, -2
        elif aa1 == aa2:
            score += 1 # match, +1
        else:
            score -= 3 # mismatch, -3

    return (score, aligned1, aligned2)

def initialize_matrix(width, height):
    """
    Description:
        This funtion initializes the alignment matrix with the first row and
        column set with the gap penalties.

    Parameters:
        - width (int): the width of the matrix
        - height (int): the height of the matrix

    Returns:
        - tbl list(list): 2D matrix to be filled in with alignment scores
    """

    # initialize table with all None values
    tbl = [[None for i in range(width + 1)] for x in range(height + 1)]
    tbl[0] = [-2 * x for x in range(width + 1)] # make first row gap penalties

    # make first column gap penalties
    counter = 0
    for row in tbl:
        row[0] = counter * -2
        counter += 1

    return tbl

def open_file(file_path):
    """
    Description:
        This function opens a protein sequence file and returns the protein
        sequence as a string.

    Parameters:
        - file_path (str): file path of the sequence file.

    Returns:
        - protein_str (str): the protein sequence from the file.
    """
    # regex to filter out non-sequence related text
    bases = re.compile("[^ATCGatcg]")

    protein_str = ""

    # open file, and extract the sequence data
    with open(file_path) as file:
        li = file.readline()

        while li:
            skip_line = bases.match(li)

            if skip_line == None:
                # skip_line == None if line only contains bases
                protein_str = li.replace("\n", "")
                return protein_str
            else:
                # was a comment or non-valid sequence line
                li = file.readline()

    return protein_str

def open_matches_file(matches_file):
    """
    Description:
        This function extracts the anchored alignment information from a matches
        file.

    Parameters:
        - matches_file (str): file path of the matches file containing anchored
            alignment information.

    Returns tuple(cov1_to_int, cov2_to_int):
        - cov1_to_int (list(list)): list of start and end positions of cov1
            anchored alignments
        - cov2_to_int (list(list)): list of start and end positions of cov2
            anchored alignments
    """
    cov1 = []
    cov2 = []

    # open file, and add build sequence data
    with open(matches_file) as file:
        li = file.readline() # skip the header line
        li = file.readline()
        while li:
            li = li.replace("\n","")
            split_anchors = li.split("\t")
            cov1.append([split_anchors[0], split_anchors[1]]) # grab cov1 data
            cov2.append([split_anchors[2], split_anchors[3]]) # grab cov2 data
            li = file.readline()

    # convert string data to ints
    cov1_to_int = [[int(pair[0]),int(pair[1])] for pair in cov1]
    cov2_to_int = [[int(pair[0]),int(pair[1])] for pair in cov2]

    return (cov1_to_int, cov2_to_int)

def save_file(score, aligned1, aligned2, seq1_file, seq2_file):
    """
    Description:
        This function is responsbile for saving the score, and aligned sequences
        to file (output.txt).

    Parameters:
        - score (int): the alignment score of aligned1 and aligned2
        - aligned1 (str): seq1 adjusted to be aligned with seq2
        - aligned2 (str): seq2 adjusted to be aligned with seq1
        - seq1_file (str): file path of first sequence file
        - seq2_file (str): file path of second sequence file
    """

    with open('output.txt', 'w') as file:
        file.write(f'Alignment score: {score}\n')
        file.write(f'{seq1_file} aligned sequence: {aligned1}\n') # write aligned1
        file.write(f'{seq2_file} aligned sequence: {aligned2}\n') # write aligned2

if __name__ == '__main__':
    # extract command line args for seq1 file, seq2 file, and [matches] file
    try:
        seq1_file = sys.argv[1] # grab sequence 1
        seq2_file = sys.argv[2] # grab sequence 2

        if len(sys.argv) > 3:
            # matches file included
            main(seq1_file, seq2_file, sys.argv[3])
        else:
            # matches file not included
            main(seq1_file, seq2_file)

    # any exceptions occured. Report error message.
    except Exception as e:
        print(f'align_proteins.py quit due to the following error: {e}')
