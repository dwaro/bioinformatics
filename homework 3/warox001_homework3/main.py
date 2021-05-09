import sys
import distances
import buildTree
import bootstrap

def main(file_name):
    '''
    Description:
        The is the main function for orchestrating homework 3. This function makes sure the data is opened,
        the original distance matrix is created, the original phylogenic tree is created, and the boostrap
        inferences of the tree are performed.
    
    Parameters:
        - file_name (str): the path to the .fna file of sequence tips for our tree (hw3.fna)
    '''

    # Question 1
    # ----------
    # open the hw3.fna file,
    # compute the distances between sequences,
    # & save the distances to genetic-distances.txt
    data = distances.open_file(file_name)
    distance_matrix = distances.compute_distances(data, True)

    # Questions 2 & 3
    # ---------------
    # creates the edges.txt file and the newick tree file (tree.txt)
    tree = buildTree.resolve_tree(distance_matrix.copy(), [str(row[0]) for row in data], True)

    # Question 4
    # ----------
    # tree image built with R script provided, saved as question_4_tree.png

    # Question 5
    # ----------
    # perform the bootstrap iterations on the tree
    # creates boostrap.txt, the corresponding boostrap confidences
    bootstrap.generate_bootstraps(tree, data)
    
if __name__ == '__main__':

    # did not pass correct number of command line args
    if (len(sys.argv) < 2):
        print('Quitting... You must provide an .fna file')

    # run program
    else:
        main(sys.argv[1])
