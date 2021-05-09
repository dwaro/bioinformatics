from align_proteins import *
import random

def main(seq1_file, seq2_file):
    """
    Description:
        This function randomly shuffles seq1_file, and then performs needleman_wunsch
        to align seq1_file and seq2_file.

    Parameters:
        - seq1_file (str): this is the file path to sequence 1
        - seq2_file (str): this is the file path to sequence 2
    """
    seq1 = open_file(seq1_file)
    seq2 = open_file(seq2_file)

    with open('permutations_S_proteins.csv', 'w') as file:
        file.write('score\n')
        for i in range(100):
            if i % 10 == 0 and i != 0: print(f'{i}%') # track progress

            # shuffle the first sequence
            new_seq_1 = list(seq1)
            random.shuffle(new_seq_1)
            new_seq_1 = ''.join(new_seq_1)

            # align sequences and write the score to file
            score, aligned1, aligned2 = needleman_wunsch(new_seq_1, seq2)
            file.write(f'{score}\n')

if __name__ == '__main__':
    # extract command line args for seq1 file, seq2 file, and [matches] file
    try:
        seq1_file = sys.argv[1]
        seq2_file = sys.argv[2]

        main(seq1_file, seq2_file)

    # issues finding command line arguments
    except Exception as e:
        print(f'permute.py quit due to the following error: {e}')
