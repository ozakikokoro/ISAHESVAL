# This python script assesses in intermediate4.txt by
# spliceAI score (acceptor increase, decrease, donor increase, decrease)
# its coordinate changers (like (+)1, 0, -43, -3).
# compare the four coordinates (only for the site with 0.5 or higher) and the original variant coordinate
# to make a region (1-based) start and end by taking the minimum and maximal coordinate.
# Finally this script add the two columns start(14th column) end (15th column) to the original imput file(13 columns)

import sys

# Get the input and output file paths from the command-line arguments
input_file_path_compare = sys.argv[1]
output_file_path_compare = sys.argv[2]

# open the input and output files
with open(input_file_path_compare, 'r') as input_file, open(output_file_path_compare, 'w') as output_file:
    # process each line of the input file
    for line in input_file:
        # split the line into columns
        columns = line.strip().split()

        # extract the relevant columns as floats
        col6, col7, col8, col9, col10, col11, col12, col13, col2 = map(float, columns[5:13] + [columns[1]])

        # create tmp1, tmp2, tmp3, and tmp4 based on the column values
        tmp1 = col2
        tmp2 = col2
        tmp3 = col2
        tmp4 = col2

        if col6 >= 0.5:
            tmp1 = col10 + col2

        if col7 >= 0.5:
            tmp2 = col11 + col2

        if col8 >= 0.5:
            tmp3 = col12 + col2

        if col9 >= 0.5:
            tmp4 = col13 + col2

        # choose the smallest value among tmp1, tmp2, tmp3, tmp4, and col2
        min_val = int(min(tmp1, tmp2, tmp3, tmp4, col2))

        # choose the largest value among tmp1, tmp2, tmp3, tmp4, and col2
        max_val = int(max(tmp1, tmp2, tmp3, tmp4, col2))

        # write the output line to the output file
        output_file.write(line.strip() + f' {min_val} {max_val}\n')
