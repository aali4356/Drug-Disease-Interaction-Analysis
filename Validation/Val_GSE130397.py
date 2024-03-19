from geode import *
import numpy as np
import os
from pprint import pprint

# Define the path to your file
# Update the path to where your combined_data_rounded.txt is located
file_path = 'GSE130397_Processed.txt'

# Initialize lists to hold gene names and the matrix data
genes = []
mat = []

# Open and read the file
with open(file_path, 'r') as f:
    # Read the first line to get the column labels
    col_labels = ['2' if x == '1' else '1' for x in f.readline().strip().split('\t')[0:]]  # Adjust the labels    
    # Read the rest of the lines for expression data
    for line in f:
        sl = line.strip().split('\t')
        gene = sl[0]
        row = list(map(float, sl[1:]))  # Convert expression values to floats
        genes.append(gene)
        mat.append(row)

# Convert the list to a numpy array
mat = np.array(mat)
print(col_labels)
print(genes)
# Since col_labels are read from the file, ensure they are correctly defined here
# This example assumes the first row in your file correctly indicates group membership
# If your groups are not simply '1' and '2', adjust col_labels accordingly

## Compute characteristic direction
chdir_res = chdir(mat, col_labels, genes, calculate_sig=0, nnull=100)
# pprint(chdir_res)
print(chdir_res[:10])
filtered_chdir_res = [(value, gene) for value, gene in chdir_res if abs(value) < 0.001]
# Filter based on a condition, e.g., abs(value) < 0.01
temp_filtered_chdir_res = [(value, gene) for value, gene in chdir_res if abs(value) < 0.001]

# Now, ensure you only keep the last 10000 entries from this filtered list
# filtered_chdir_res = temp_filtered_chdir_res[-10000:]

print(filtered_chdir_res[:10])
# Saving results
# Define the path for the output text file
output_file_path = 'chdir_GSE130397_results.txt'

# Open the output file in write mode
with open(output_file_path, 'w') as out_file:
    # Iterate over each gene-value pair in chdir_res
    for value, gene in chdir_res:
        # Write the gene and value to the file, formatted as "Gene, value"
        out_file.write(f"{gene}, {value}\n")

print(f"Results saved to {output_file_path}")
