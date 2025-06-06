# This python script make intermediate2.txt into chr, coordinate, ref, alt (splitted),Gene_symbol,
# spliceAI score (acceptor increase, decrease, donor increase, decrease) and its corresponding positions
# (as above (position to the coordinate)) space delmited file (13 columns).

import sys

# Get the input and output file paths from the command-line arguments
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

# Open the input and output files
with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    # Loop through each line in the input file
    for line in input_file:
        # Split the line into columns using space as a delimiter
        columns = line.strip().split()
        # Separate the 5th column using "|" as a delimiter
        fifth_column = columns[4].split('|')
        # Combine the first three columns and the separated 5th column into a single list
        new_columns = columns[:3] + fifth_column
        # Write the modified line to the output file, using space as a delimiter between columns
        output_file.write(' '.join(new_columns) + '\n')
