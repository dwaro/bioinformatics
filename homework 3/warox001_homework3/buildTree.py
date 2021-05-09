import copy, collections
from operator import itemgetter

def resolve_tree(distance_matrix, sequence_list, save):
    '''
    Description:
        This function is responsible for generating a tree given a distance matrix and a list of sequence tips.
        This function performs this operation by performing a neighbor joining algorithm

    Parameters:
        - distance_matrix (list(list)): matrix of distances between all tip pairs
        - sequence_list (list): list of all the sequence tips in the data
        - save (boolean): true/false value that decides if we want to save the edges & newick tree files.

    Returns:
        - tree (dict): the representative newick tree given the distance matrix. 
    '''
    tree = {}
    size = len(sequence_list) # will mutate as we progress
    sizeCopy = size # want to keep this unchanged
    count = size * 2 - 2 # number of initial tips + number of internal nodes that will be formed
    lookup = {sequence_list[i]: i + 1 for i in range(size)} # mapping of each tip id to it's position in the sequence list

    # Neighbor joining algorithm, completes size - 2 joins
    while count > sizeCopy:

        # create the space for the q-matrix
        q_matrix = [[None for seq in sequence_list] for sequence in sequence_list]

        min_val = size - 1 # track the minimum val
        min_row = 0 # track the row of the minimum val
        min_col = 0 # track the column of the minimum val

        # calculate q-matrix
        # we only need to do 1/2 the matrix by using the fact that the matrix is symmetrical
        # across the diagonal
        for row in range(size):
            for col in range(size - row):
                # skip self comparison
                if row == col + row:
                    val = None
                else:
                    # calculate distance
                    val = (size - 2) * distance_matrix[row][col+row] - sum(distance_matrix[row]) - sum([sequence[col+row] for sequence in distance_matrix])
                    
                    # update minimum value accordingly
                    if (val < min_val):
                        min_val, min_row, min_col = (val, row, col + row)

                # update q-matrix value
                q_matrix[row][col+row] = val
                q_matrix[col+row][row] = val

        # count is the number of the new node
        # calculate distances from the row tip/node to the new node, and the distance from the col tip/node to the new node
        dist_row_u = 0.5 * distance_matrix[min_row][min_col] + (1 / (2 * (size - 2))) * (sum(distance_matrix[min_row]) - sum([seq[min_col] for seq in distance_matrix]))
        dist_col_u = distance_matrix[min_row][min_col] - dist_row_u

        # grab the node position value, and then update the tree with the two descendants
        actual_row = sequence_list[min_row] if 62 <= int(sequence_list[min_row]) <= 120 else lookup[sequence_list[min_row]]
        actual_col = sequence_list[min_col] if 62 <= int(sequence_list[min_col]) <= 120 else lookup[sequence_list[min_col]]
        tree[count] = [(count, actual_row, dist_row_u), (count, actual_col, dist_col_u)]

        # get the distances to the new node from the rest of the tree
        new_node_distances = [0.5 * (distance_matrix[min_row][i] + distance_matrix[i][min_col] - distance_matrix[min_row][min_col]) for i in range(size) if i not in (min_row, min_col)]
        new_node_distances.insert(0, 0)

        # remove the old nodes from the table, and insert the newly formed node
        sequence_list.pop(min(min_row, min_col))
        sequence_list.pop(max(min_row, min_col)-1)
        sequence_list.insert(0, count)

        # update the distance matrix, and add the new node distances to the front of the matrix
        distance_matrix = [distance_matrix[row] for row in range(size) if row not in (min_row, min_col)]
        for row in distance_matrix:
            row.pop(min(min_row, min_col))
            row.pop(max(min_row, min_col)-1)

        # update the matrix
        distance_matrix.insert(0, new_node_distances)
        counter = 1
        for row in distance_matrix[1:]:
            row.insert(0, new_node_distances[counter])
            counter += 1

        # update trackers
        size -= 1
        count -= 1

    # add the final connection in the tree
    lst = tree[sizeCopy+1]
    # if the remaining connection is to a newly created internal node (int)
    if type(sequence_list[len(sequence_list)-1]) == int:
        lst.append((sizeCopy+1, sequence_list[len(sequence_list)-1], distance_matrix[0][1]))

    # else last connection is to an original tip (str)
    else:
        lst.append((sizeCopy+1, lookup[sequence_list[len(sequence_list)-1]], distance_matrix[0][1]))
    tree[sizeCopy+1] = lst

    # sort the children at each node on the descendant node value for consistency
    try:
        tree = {node: sorted(tree[node], key=itemgetter(1)) for node in tree}
    except:
        print(tree)

    # save the edge and newick tree files
    if save: save_files(tree, sizeCopy, lookup)

    return tree

def save_files(tree, size, lookup):
    '''
    Description:
        This function is responsible for facilitating the saving of the edges and the newick tree files.

    Parameters:
        - tree (dict): the dictionary representing the tree. node: [(node, descendant1, distance1), (node, descendant2, distance2)]
        - size (int): the number of tips
        - lookup (dict): a mapping of sequence_tip_id: position in data (1 - 61)
    '''
    root_children = tree[size+1] # start at the root
    root_children.reverse() # trying to match solution's order

    # save the edges file
    with open('edges.txt', 'w') as file:
        # save each branch from root
        for child in root_children:
            save_edges(tree, child, file, size)

    # save the newick tree
    position_map = {lookup[key]: key for key in lookup} # inverse of lookup, tip_position: tip_id
    final_str = ''
    for child in root_children:
        final_str += save_newick(tree, child, position_map) + ','
    final_str = f'({final_str[:len(final_str) - 1]});'

    with open('tree.txt', 'w') as file:
        file.write(final_str)
    
def save_edges(tree, node, file, start):
    '''
    Description:
        This funciton writes the tree's edges to file using a recursive approach.

    Parameters:
        - tree (dict): the dictionary representing the tree. node: [(node, descendant1, distance1), (node, descendant2, distance2)]
        - node (int): the current node to continue recursing through the tree's branches
        - file (file): the file currently being written to
        - start (int): the node to treat as the root in the tree (62)
    '''
    # preorder, first write to file
    file.write(str(node[0]) + '\t' + str(node[1]) + '\t' + str(node[2]) + '\n')

    # recurse if the node is an internal node (original tips won't have children)
    if node[1] > start:
        # recurse for each child (left and right)
        for child in tree[node[1]]:
            save_edges(tree, child, file, start)

def save_newick(tree, ancestor, pos_map):
    '''
    Description:
        This function builds a branch from root and adds it to the final string to be saved representing
        the newick tree. (Recursive)

    Parameters:
        - tree (dict): the dictionary representing the tree. node: [(node, descendant1, distance1), (node, descendant2, distance2)]
        - ancestor (int): the previous ancestor node
        - pos_map (dict): a mapping of sequence_tip_id: position in data (1 - 61)
    
    Returns:
        - (str) (lef subtree, right subtree):ancestor_distance
    '''
    current = ancestor[1] # current node in tree

    # stop condition, we've reached a tip, no longer recurse as tips don't have children
    if current not in tree:
        return f'{pos_map[current]}:{ancestor[2]}'

    # recurse deeper, moving down the left and right subtrees
    children = tree[current]
    left = save_newick(tree, children[1], pos_map)
    right = save_newick(tree, children[0], pos_map)
    
    # postorder
    return f'({left},{right}):{ancestor[2]}'