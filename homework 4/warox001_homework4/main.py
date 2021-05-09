import sys
import matplotlib.pyplot as plt
from operator import itemgetter

def main(sequence_file):
    '''
    Description:
        Facilitates the execution of calculating/plotting conservation rates and
        calculating/plotting variable regions.
    
    Parameters:
        - sequence_file (str): the file path of the sequence-with-primers.fna file containing aligned sequences
    '''
    # Steps 1 & 2
    conservation_rates = calculate_conservation_rate(sequence_file) # open data, and calculate conservation rates
    smoothed_data = smooth_data(conservation_rates, 6) # smooth the data
    plot_conservation_rates(smoothed_data, len(conservation_rates[0])) # plot just the conservation rates

    # Steps 3 & 4
    variable_regions = generate_variable_regions(conservation_rates, 30, 100, num_regions=9) # calculate the variable regions
    plot_conservation_rates(smoothed_data, len(conservation_rates[0]), variable_regions=variable_regions) # plot var regions over conservation rates

def calculate_conservation_rate(in_file):
    '''
    Description:
        This function opens the sequence-with-primers data and calculates the conservation rates.
    
    Parameters:
        - in_file (str): the file path of the sequence-with-primers.fna file containing aligned sequences
    
    Returns:
        - output ([list(),list()]): list(sequence_positions, conservation_rates)
    '''
    seq_size = 0

    # get the length of sequences
    with open(in_file, 'r') as file:
        line = file.readline()

        while '>' in line:
            line = file.readline()

        line = line.replace('\n', '')
        seq_size = len(line)

    # data to hold the counts of each ACTG base at each position
    data = {num: {'A':0,'C':0,'T':0,'G':0} for num in range(1, seq_size + 1)}

    # this will be the number of sequences in the file
    num_sequences = 0

    # open the file and count ACTG occurences at each position of each sequence
    with open(in_file, 'r') as file:
        line = file.readline()

        while line:
            # check if it's the sequence line
            if '>' not in line:
                num_sequences += 1

                # record the ACTG base at each position, ignore anything besides ACTG
                for i in range(seq_size):
                    if line[i] in ('A', 'C', 'T', 'G'):
                        # update base count
                        current = data[i+1]
                        current[line[i]] = current[line[i]] + 1
                        data[i+1] = current
                
            line = file.readline() # move forward in file
    
    output = [[i for i in range(1, seq_size+1)],[]]
    rates = []

    # write the conservation rates to 'solution-problem-1.txt' file
    with open('solution-problem-1.txt', 'w') as out_file:
        keys = sorted(list(data.keys())) # increasing order of positions for saving to file

        # for each position, calculate the conservation rate
        for key in keys:
            dic = data[key] # position_number == key: {A: a_count, C: c_count, T: t_count, G: g_count }
            conservation_rate = max(dic.values()) / num_sequences # calc conservation rate
            rates.append(conservation_rate) # record the rate
            out_file.write(str(conservation_rate) + '\n') # write to file
    
    output[1] = rates
    return output

def smooth_data(conservation_rates, step):
    '''
    Description:
        This function smooths the data, to make it more simplified for plotting purposes.
    
    Parameters:
        - conservation_rates ([list(), list()]): the conservation rates at each position
        - step (int): the step size in between points plotted, and the window size will be 2 * step + 1

    Returns:
        - [smoothed_x (positions), smoothed_rates]: the conservation rates smoothed using a sliding window average
    '''
    raw_rates = conservation_rates[1]

    smoothed_x = [] # the positions that will actually get plotted
    smoothed_rates = [] # the average rate for a window size: step * 2 + 1 (+/- step size)

    # calculate smoothed rates
    for i in range(step, len(conservation_rates[0]) - step, step):
        smoothed_x.append(i+1) # record the position
        window_avg = sum(raw_rates[i-step:i+step+1]) / (step * 2 + 1) # sliding window average
        smoothed_rates.append(window_avg)

    return [smoothed_x, smoothed_rates]

def plot_conservation_rates(conservation_rates, size, variable_regions=None):
    '''
    Description:
        This function plots the conservation rates & variable regions
    
    Parameters:
        - conservation_rates ([list(), list()]): the conservation rates at each position
        - size (int): the number of positions (columns) in the sequences
        - variable_regions (optional, default: None): the data on the variable_regions (list([region_variance, region_start, region_end]))
    '''
    # plot(x series, y series)
    # ------------------------
    plt.figure(figsize=(15,8)) # figure size
    plt.plot(conservation_rates[0], conservation_rates[1], color='red', linewidth=2) # plot smoothed conservation rates

    # set axes & labels & grid
    plt.xlim([0, size])
    plt.ylim([0, 1])
    plt.ylabel('% Sequence Identity', weight='bold', size='large')
    plt.xlabel('Gene Position', weight='bold', size='large')
    plt.grid(linestyle="--", linewidth=1)

    # save figure (problem 2)
    plt.savefig('solution-problem-2.pdf', bbox_inches='tight')

    # if also plotting variable regions...
    if variable_regions is not None:

        # plot each variable region
        for region in variable_regions:
            x_series = [] # the gene positions

            # grab all the gene positions for this region
            start, end = region[1], region[2]
            for i in range(start, end+1):
                x_series.append(i)

            y_constant = [0.5 for value in x_series] # plot everything to y-value 0.5
            plt.plot(x_series, y_constant, color='black', linewidth=4)
        
        # save figure (problem 4)
        plt.savefig('solution-problem-4.pdf', bbox_inches='tight')

def generate_variable_regions(conservation_rates, min_region, max_region, num_regions=9):
    '''
    Description:
        This function calculates 9 variable regions. The method it employs is calculating variances for all
        sequential subarray gene positions for the entire 1514 positions. There are two requirements for
        calculating the variances of a subarray, it must be size min_region or larger, and have a conservation
        rate range of at least 0.2 in its region.
    
    Parameters:
        - conservation_rates ([list(), list()]): the conservation rates at each position
        - min_region (int): the minimum size/#_of_gene_positions to use for a variable region. Also used as 
            the threshold for merging close regions.
        - max_region (int): the maximum region size we want to consider.

    Returns:
        - var_regions (list([region_variance, region_start, region_end], ... , region 9)): list of the variable regions
    '''
    # all of the sequential gene position subarray regions of min_region size or greater
    sub_var_regions = []

    # raw rates data
    og_data = conservation_rates[1]
    size = len(conservation_rates[0]) # number of gene positions

    # calculate variances for all sequential slices of the conservation rates, with a min. region slice of min_region
    for i in range(len(conservation_rates[0])):
        # making an assumption that a variable region needs to be at least min_region bases long
        tmp_size = min_region

        # for all sub-array sizes possible between i & the size (1514 gene positions)
        while (tmp_size < max_region and i + tmp_size <= size):
            data = og_data[i:i+tmp_size] # get the subarray (slice) of the gene positions
            mean = sum(data) / len(data) # calculate the mean of the region
            variance = sum([(val - mean)**2 for val in data]) / (len(data) - 1) # calc variance

            # try to enforce some level of variation by requring conservation rate range > 0.2
            if max(data) - min(data) > 0.25:
                sub_var_regions.append([variance, i, i + tmp_size - 1])
            tmp_size += 1

    # sort the variable regions by variance decreasing
    sub_var_regions.sort(key=itemgetter(0), reverse=True)
    
    # find the 9 most variable, distinct, regions
    var_regions = [sub_var_regions[0]] # init with most variable region
    while len(var_regions) < num_regions:
        position = 1 # tracker

        # get the 9 most variable regions that don't overlap
        while len(var_regions) < num_regions and position < len(sub_var_regions):
            region = sub_var_regions[position]
            start, end = region[1], region[2]
            valid = True

            # make sure there is no overlap
            for variable_region in var_regions:
                if end < variable_region[1] or start > variable_region[2]:
                    # no overlap
                    pass
                else:
                    # overlap
                    valid = False
                    break
            
            # if a valid region, add it, and remove it from sub_var_regions
            if valid:
                var_regions.append(region)
                sub_var_regions.pop(position)
                position -= 1 # account for removal from sub_var_regions

            position += 1

        # if any of the 9 regions are close together, merge them together
        # and go find the next most variable region
        merging = True
        while merging:
            merging = False

            for i in range(len(var_regions)):
                start, end = var_regions[i][1], var_regions[i][2]
                
                # compare this.region to the subsequent regions
                for j in range(i+1,len(var_regions)):
                    # if this.region is close to any of the other regions, merge them
                    if start > var_regions[j][2] and start - var_regions[j][2] < 20:
                        merging = True
                        var_regions.append([None, var_regions[j][1], end]) # merged region
                        var_regions.pop(i) # remove this.region
                        var_regions.pop(j-1) # remove the region it was joined with
                        break
                    elif end < var_regions[j][1] and var_regions[j][1] - end < 20:
                        merging = True
                        var_regions.append([None, start, var_regions[j][2]]) # merged region
                        var_regions.pop(i) # remove this.region
                        var_regions.pop(j-1) # remove the region it was joined with
                        break
                
                if merging:
                    break
        
    # save the variable regions to 'solution-problem-3.txt'
    with open('solution-problem-3.txt', 'w') as out_file:
        for region in var_regions:
            out_file.write(f'{(region[1]+1)}\t{(region[2]+1)}\n') # adjust for 1-based indexing

    return var_regions

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print('You must provide 1 argument (the name of the sequence file) to the program')
    
    else:
        main(sys.argv[1])