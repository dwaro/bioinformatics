import random
import distances
import collections
import buildTree

def generate_bootstraps(tree, data):
    '''
    Description:
        This function performs 100 inferences of bootstrap sampling to generate a partition
        confidence value for each node in the original phylogenic tree that we constructed.

    Parameters:
        - tree (dict): the original pyhlogenic tree constructed for our data
        - data (list(list)): the sequence tips data in the format [[tip_id, sequence], ... , [last_tip_id, last_sequence]]
    '''
    length = len(data[0][1]) # length of a sequence (they're all the same length)

    # this will hold the count of how often a node's substructure is found in the boostrapped trees
    node_counts = {node:0 for node in tree}

    # this is a list of all the possible sub-tree tips from the original tree
    # keeping it sorted for easy comparison later
    possible_partitions = {}
    for node in tree:
        node_level = []
        for child in tree[node]:
            tip_list = []
            get_tips(tree, child[1], tip_list, len(data) + 1)
            node_level.append(sorted(tip_list))
        possible_partitions[node] = sorted(node_level)
    
    # 100 inferences of the tree
    for inference in range(100):
        # temporary storage of the sequences for each node
        temp = {seq[0]:'' for seq in data}

        # sample length columns
        for step in range(length):
            # the random sequence position to be added to the bootstrap sample
            col = random.randint(0, length - 1)
            
            # add the base at position col of each original sequence to the boostrap sequence for each tip
            for sequence in data:
                val = sequence[1][col] # grab the base at position: col
                node = sequence[0] # get the node 
                temp[node] = temp[node] + val # update the bootstrap sequence

        # generate new tree for this bootstrap sample
        new_data = [[key, temp[key]] for key in temp]
        new_distance_matrix = distances.compute_distances(new_data, False)
        new_tree = buildTree.resolve_tree(new_distance_matrix, [str(row[0]) for row in new_data], False)

        # compare the new tree to the old tree and count bootstrap confidence
        for node in new_tree:
            node_level = []
            for child in new_tree[node]:
                tip_list = []
                get_tips(tree, child[1], tip_list, len(data) + 1)
                node_level.append(sorted(tip_list))
            node_level = sorted(node_level)

            # compare new_tree's tip partitions for node to all the possible partition from original tree
            for og_node in possible_partitions:
                if possible_partitions[og_node] == node_level:
                    node_counts[og_node] = node_counts[og_node] + 1 # found a match

    # adjust the counts to be 0 - 1
    node_counts = {node:(node_counts[node] / 100) for node in node_counts}
    
    # save the bootstrap counts for each node
    with open('bootstrap.txt', 'w') as file:

        visited = []

        # save in order as they appear in the edges.txt file
        with open('edges.txt', 'r') as edges:

            line = edges.readline()
            while line:
                node = int(line.split('\t')[0])

                # if it's an internal node
                if node > len(data):
                    # if we haven't visited this node yet
                    if node not in visited:
                        val = node_counts[int(node)] # get the count
                        if val == 0.0: val = 0
                        if val == 1.0: val = 1
                        file.write(f'{val}\n')
                        visited.append(node)

                line = edges.readline()

def get_tips(tree, current, tip_list, size):
    '''
    Description:
        This function generates the tip partitions recursively for a given node

    Parameters:
        - tree (dict): the original tree created
        - current (int): the current node we're at during our search down the tree
        - tip_list (list): the list of tips for a given start node
        - size (int): the number of tips + 1 from the data
    '''
    # check if it's a tip node, if so add the tip_list
    if current < size:
        tip_list.append(current)
    
    # if it's not a tip node, continue recursion looking for tips
    else:
        for child in tree[current]:
            get_tips(tree, child[1], tip_list, size)