from geode import *
import numpy as np
import os
from pprint import pprint

# Define the path to your folder containing the input files
folder_path = 'Data/Processed'

# Define the path to the output directory
output_dir = 'Data/Data_Output'

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Iterate over each file in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('.txt'):  # Process only .txt files
        file_path = os.path.join(folder_path, file_name)
        
        # Remove the "Processed.txt" part from the file name
        output_file_name = file_name.replace('_Processed.txt', '')
        
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
        # print(mat)
        # Convert the list to a numpy array
        mat = np.array(mat)
        
        print(col_labels)
        print(genes)
        
        # Compute characteristic direction
        chdir_res = chdir(mat, col_labels, genes, calculate_sig=0, nnull=100)
        
        print(chdir_res[:10])
        
        filtered_chdir_res = [(value, gene) for value, gene in chdir_res if abs(value) < 0.001]
        
        print(filtered_chdir_res[:10])
        
        # Saving results
        # Define the path for the output text file
        output_file_path = os.path.join(output_dir, f'chdir_{output_file_name}_results.txt')
        
        # Open the output file in write mode
        with open(output_file_path, 'w') as out_file:
            # Iterate over each gene-value pair in chdir_res
            for value, gene in chdir_res:
                # Write the gene and value to the file, formatted as "Gene, value"
                out_file.write(f"{gene}, {value}\n")
        
        print(f"Results for {output_file_name} saved to {output_file_path}")