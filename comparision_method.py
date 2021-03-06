#!/usr/bin/env python

################################################################################
# Compute the pearson correlation between two base-pairing
# probibility matricies from RNA structure 
#
# Author - Kobie Kirven 
################################################################################


# Imports
import argparse
import numpy as np 
import math
import scipy.stats as st

# Add the command-line arguments
argparser = argparse.ArgumentParser(description='Compute the pearson correlation betwen two base-pairing probibility matricies from RNA structure')
argparser.add_argument('-1', '--file1', help='BP File 1', 
dest="file1", required=True)
argparser.add_argument('-2', '--file2', help='BP File 2', dest="file2", required=True)
argparser.add_argument('-c', '--comp', help="Comparison method: pearson, pearson-not-sum, rmsd", dest="comp", required=True)
args = argparser.parse_args()

# Functions
#----- 
# pearson-sum 
# pearson-not-sum
# rmsd 

def turn_propibility_plot_to_matrix(fileName):
    """
    Turn the probibility of base pairing into
    a numpy matrix 

    parameters: 
     - fileName - the name of the file to read from
    """
    # Read the file
    with open(fileName) as fn:
        lines = fn.readlines()
        matrix = np.zeros((int(lines[0].strip("\n")), int(lines[0].strip("\n"))))
        for line in lines[2:]:
            line = line.strip("\n").split("\t")
            matrix[int(line[0])-1, int(line[1])-1] = math.pow(10, -1 * float(line[2]))
        return matrix 

def get_proibilites(fileName):
    """
    Get the probibility of base pairing from a file
    """
    #
    # with open(fileName) as fn:
    fn = open(fileName)
    lines = fn.readlines()
    base_pairs = []
    probs = []

    for line in lines[2:]:
        line = line.strip("\n").split("\t")
        base_pairs.append([int(line[0]), int(line[1])]) 
        probs.append(math.pow(10, -1 * float(line[2])))

    return base_pairs, probs

def compute_pearson_corr(matrix1, matrix2):
    """
    Compute the pearson correlation between two bpp vectors
    """
    vector1 = np.sum(matrix1, axis=1)
    vector2 = np.sum(matrix2, axis=1)

    # Compute the pearson correlation of the bpp vectors
    return st.pearsonr(vector1, vector2)[0]

def pearson_not_sum(probs1, probs2):
    """
    Compute the pearson correlation between two bpp vectors
    """
    # Compute the pearson correlation of the bpp vectors
    return st.pearsonr(probs1, probs2)[0]

def compute_rmsd(pr1, pr2):
    """
    Compute the rmsd between the two bpp matricies
    """
    diffs = []
    # Compute the rmsd of the bpp matricies
    for i in range(len(pr1)):
        diffs.append((pr1[i] - pr2[i])**2)
    return np.sqrt(sum(diffs)/len(pr1))
   


# Run the tool
if __name__ == "__main__":
    matrix1 = turn_propibility_plot_to_matrix(args.file1)
    matrix2 = turn_propibility_plot_to_matrix(args.file2)

    probs1 = get_proibilites(args.file1)
    probs2 = get_proibilites(args.file2)

    # Compute the pearson correlation
    if args.comp == "pearson":
        print(compute_pearson_corr(matrix1, matrix2))
    
    if args.comp == "pearson-not-sum":
        print(pearson_not_sum(probs1[1], probs2[1]))
    
    if args.comp == "rmsd":
        print(compute_rmsd(probs1[1], probs2[1]))
    

