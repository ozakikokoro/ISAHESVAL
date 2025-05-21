# This python script calculate the minimum number
# the output variable "min_val_from_line" can be captured in the shell script.

import sys
input = sys.argv[1]

min_val = float('inf')
min_line = ''

with open(input, 'r') as f:
    for line in f:
        cols = line.strip().split()
        val = float(cols[2])
        if val < min_val:
            min_val = val
            min_line = line.strip()

# `min_val` contains the minimum value from the third column, and `min_line`
# contains the line it was found in.

# Extract the value from `min_line` and print it with 10 decimal places
min_value_str = min_line.split()[2]
min_value_float = float(min_value_str)

print(f"{min_value_float:.10e}")
# Use formatted string literals for precision (floating point expression(scientific notation) with 10 digits after the decimal point)
