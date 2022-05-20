################################################################################
# Compute the pearson correlation betwen two base-pairing
# probibility matricies from RNA structure 
#
# Author - Kobie Kirven 
################################################################################

#!/usr/bin/env python3

# Imports
import argparse
import subprocess
import tempfile
import numpy as np 
import math
import scipy.stats as st

# Add the command-line arguments
argparser = argparse.ArgumentParser(description='Compute the pearson correlation betwen two base-pairing probibility matricies from RNA structure')
argparser.add_argument('-1', '--seq-1', help='Sequence 1', dest="seq1", required=True)
argparser.add_argument('-2', '--seq-2', help='Sequence 2', dest="seq2", required=True)
argparser.add_argument('-s', '--singularity', help='Path to singularity file', required=True)
args = argparser.parse_args()

# Functions

def turn_propibility_plot_to_matrix(fileName):
    """
    Turn the probibility of base pairing into
    a numpy matrix 
    """
    # Read the file
    with open(fileName) as fn:
        lines = fn.readlines()
        matrix = np.zeros(int(lines[0].strip("\n")), int(lines[0].strip("\n")))
        for line in lines[2:]:
            line = line.strip("\n").split(" ")
            matrix[int(line[0])-1, int(line[1])-1] = math.pow(10, -1 * float(line[2]))
        return matrix 

def get_bpp_matrix(singularity_path, sequence, temp_dir):
    """
    Get the base-pairing probibility matrix for an mRNA sequence
    """
    # Create a temporary file
    temp = tempfile.NamedTemporaryFile(delete=False)

    # Write the sequence to the temporary file
    temp.write(f">test\n{sequence}\n".encode('utf-8'))
    # Close the file
    temp.close()
    # Get the base-pairing probibility matrix
    subprocess.run(['singularity', 'exec', singularity_path, "partition", temp.name, f"{temp_dir}/out.psf"])

    # convert the psf file to the matrix file
    subprocess.run(['singularity', 'exec', singularity_path, "ProbabilityPlot", "-t", f"{temp_dir}/out.psf", f"{temp_dir}/out.txt"])

    # Delete the temporary file
    subprocess.call(['rm', temp.name])

    # Return the matrix
    return turn_propibility_plot_to_matrix(f"{temp_dir}/out.txt")

def compute_pearson_corr(matrix1, matrix2):
    """
    Compute the pearson correlation between two bpp vectors
    """
    vector1 = np.sum(matrix1, axis=1)
    vector2 = np.sum(matrix2, axis=1)

    # Compute the pearson correlation of the bpp vectors
    return st.pearsonr(vector1, vector2)[0]


# Run the tool
with tempfile.TemporaryDirectory() as tmpdirname:
    print('created temporary directory', tmpdirname)
    # Get the base-pairing probibility matricies
    
    matrix1 = get_bpp_matrix(args.singularity, args.seq1, tmpdirname)
    matrix2 = get_bpp_matrix(args.singularity, args.seq2, tmpdirname)

    # Compute the pearson correlation
    print(f"Correlation:\t{compute_pearson_corr(matrix1, matrix2)}")
    

