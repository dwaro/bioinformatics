import sys, re

def main(input_path, output_path):
  """
    Main function to count_codons. Responsible for orchestrating the opening and parsing of sequence data, creating a frequency map of codons, and outputting the counts to the specified output file.

    Parameters
    ----------
    @param input_path {str}: the file location of the input sequence data
    @param output_path {str}: the location to put the output file
  """

  sequence_data = open_file(input_path)
  codon_map = count_codons(sequence_data)
  output_frequencies(codon_map, output_path)

  return "done"

def open_file(input_path):
  """
    This function is responsible for opening the input data file, and parsing the sequence data from it.

    @param input_path {str}: the file location of the input sequence data
    @returns {str}: the combined string of all the sequence data
  """

  sequence_data = []

  # regex to filter out non-sequence related text
  bases = re.compile("[^ATCGatcg]")

  # open file, and add build sequence data
  with open(input_path) as in_file:
    li = in_file.readline()
    count = 0
    while li:
      skip_line = bases.match(li)
      if skip_line == None:
        sequence_data.append(li.replace("\n", ""))
      li = in_file.readline()

  return "".join(sequence_data)

def count_codons(sequence_data):
  """
    This function creates a python dictionary / mapping of the number of each codon in the sequence data.

    @param sequence_data {str}: the sequence data string to be parsed.
    @returns codon_map {dict}: the counts of each codon.
  """

  codon_map = {}
  stop, size = 3, len(sequence_data)

  # step through sequence with a window size of 3 to count codons
  while stop <= size:
    codon = sequence_data[stop - 3 : stop] # window of 3 residues (1 codon)

    if codon not in codon_map:
      codon_map[codon] = 1
    else:
      codon_map[codon] = codon_map[codon] + 1

    stop += 3 # move window to next codon
  
  return codon_map
  
def output_frequencies(codon_map, output_path):
  """
    This function writes the counts of codon frequencies to an output .csv file in the format: codon,frequency\n

    @param codon_map {dict}: the counts of each codon.
    @param output_path {str}: the location to put the output file
  """

  # ensure .csv extension
  if '.csv' not in output_path: output_path += '.csv'

  keys = list(codon_map.keys())
  
  # get the last codon, so way can check not to add extra empty line to the end
  end = keys[len(keys)-1]

  # write file
  with open(output_path, 'w') as output:
    for codon in codon_map:
      li = codon + ',' + str(codon_map[codon])
      if codon != end: li += '\n'
      output.write(li)

if __name__ == "__main__":
  main(sys.argv[1], sys.argv[2])