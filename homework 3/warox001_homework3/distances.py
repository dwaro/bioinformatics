def open_file(file_name):
    '''
    Description:
        This function is responsible for opening the list of sequence tips data.

    Parameters:
        - file_name (str): the path to the .fna file of sequence tips for our tree (hw3.fna)

    Returns:
        - data (list(list)): the opened sequence tips data in the format [[tip_id, sequence], ... , [last_tip_id, last_sequence]]
    '''
    data = []

    with open(file_name, 'r') as file:
        tmp = []
        line = file.readline()
        while line:

            # grab the tip id
            if '>' in line:
                line = line.replace('>', '')
                line = line.replace('\n', '')
                tmp.append(line)

            # grab the sequence
            else:
                line = line.replace('\n', '')
                tmp.append(line)
                data.append(tmp.copy())
                tmp = []
            line = file.readline()
    
    return data

def compute_distances(data, save):
    '''
    Description:
        This function is responsible for computing the genetic distance matrix between all of the sequences. It then
        makes sure to save the matrix to genetic-distances.txt

    Parameters:
        - data (list(list)): the sequence tips data in the format [[tip_id, sequence], ... , [last_tip_id, last_sequence]]

    Returns:
        - matrix (list(list)): the distance_matrix holding genetic distances between each pair of sequences
    '''
    size = len(data)

    # this will serve as the genetic distance matrix
    matrix = [[None for x in range(size)] for y in range(size)]

    # iterate over all sequences and calculate the distance between each pair
    for i in range(size):
        for j in range(size):

            # calculate distance
            if i != j:
                seqLength = len(data[i][1])
                same = 0

                # calculate the genetic distance between sequence i and sequence j
                for t in range(seqLength):
                    if data[i][1][t] == data[j][1][t]:
                        same += 1
                
                # save sequence length
                matrix[i][j] = 1 - (same / seqLength)

            # save the diagonal as 0s
            else:
                matrix[i][j] = 0
    
    if save: save_distances(matrix, data, size) # save the distance matrix to file
    return matrix

def save_distances(matrix, data, size):
    '''
    Description:
        This function is responsible for saving the genetic-distances file

    Parameters:
        - matrix (list(list)): distance_matrix of sequence similarities
        - data (list(list)): list of the sequence tips
        - size (int): the number of sequences
    '''
    with open('genetic-distances.txt', 'w') as output:
        # write first header line
        line = '\t'.join([str(sequence[0]) for sequence in data])
        output.write('\t' + line + '\n')

        # write all the rows to the output file
        for row in range(size):
            line = str(data[row][0]) + '\t'
            line += '\t'.join([str(matrix[row][col]) for col in range(size)])
            output.write(line + '\n')