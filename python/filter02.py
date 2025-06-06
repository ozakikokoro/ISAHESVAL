# This python script filter intermediate3.txt by
# spliceAI score (acceptor increase, decrease, donor increase, decrease)
# at least one of the four scores equal to or greater than 0.2.

import sys

# Get the input and output file paths from the command-line arguments
input_file_path_filter = sys.argv[1]
output_file_path_filter = sys.argv[2]

import pandas as pd
# Read data.txt file into a pandas DataFrame
df = pd.read_csv(input_file_path_filter, sep=' ', header=None)

# Exclude lines where columns 6-13 are all "."
df = df.loc[~(df.iloc[:, 5:] == ".").all(axis=1)]

# Convert columns 6-9 to numeric values
df.iloc[:, 5:9] = df.iloc[:, 5:9].apply(pd.to_numeric, errors='coerce')

# Filter lines based on numerical values in columns 6-9
df = df.loc[(df.iloc[:, 5:9] >= 0.2).any(axis=1)]

# Write filtered data to a new file
df.to_csv(output_file_path_filter, sep=' ', header=None, index=False)
