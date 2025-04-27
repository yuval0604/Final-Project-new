import math
import sys
import pandas as pd
import numpy as np
import symnmf_capi as symnmf
from sklearn.metrics import silhouette_score


#Configuration
ERROR_MSG = "An Error Has Occurred"
CSV_SEPARATOR = ','
SEED = 1234
BETA = 0.5
EPSILON = 1e-4
MAX_ITERS = 300

np.random.seed(SEED)
#Utility Functions
def exit_with_error(msg=ERROR_MSG):
    print(msg)
    sys.exit(1)

def display_matrix(matrix):
    """Prints a matrix with values formatted to 4 decimal digits."""
    for row in matrix:
        print(','.join(f"{val:.4f}" for val in row))

def read_data_file(filepath):
    """Reads a text file into a representing a matrix."""
    if not filepath.endswith(".txt"):
        exit_with_error()
    try:
        df = pd.read_csv(filepath, delimiter=CSV_SEPARATOR, header=None)
        return df.values.tolist()
    except Exception:
        exit_with_error()

def initialize_H_matrix(laplacian_matrix, k):
    """
    Generates an initial guess for H matrix using random uniform values.

    @param laplacian_matrix: Laplacian matrix (2D list)
    @param k: Number of clusters
    @return: H matrix (2D list)
    """
    mean_val = np.mean(laplacian_matrix)
    bound = 2 * math.sqrt(mean_val / k)
    return [[np.random.uniform(0, bound) for _ in range(k)] for _ in range(len(laplacian_matrix))]

#Goals
def sym(data):
    """Computes and prints the similarity matrix."""
    similarity_matrix = symnmf.sym(data)
    display_matrix(similarity_matrix)

def ddg(data):
    """Computes and prints the diagonal degree matrix."""
    ddg_matrix = symnmf.ddg(data)
    display_matrix(ddg_matrix)

def norm(data):
    """Computes and prints the normalized Laplacian matrix."""
    norm_matrix = symnmf.norm(data)
    display_matrix(norm_matrix)

def symnmf_algo(data, k, beta, eps, max_iters):
    """
    Executes the SymNMF algorithm and prints the resulting H matrix.

    @param data: Data matrix
    @param k: Number of clusters
    """
    norm_laplacian = symnmf.norm(data)
    initial_H = initialize_H_matrix(norm_laplacian, k)
    final_H = symnmf.symnmf(initial_H, norm_laplacian, beta, eps, max_iters)
    return final_H


def main():
    if len(sys.argv) != 4:
        exit_with_error()

    try:
        k = int(sys.argv[1])
        goal = sys.argv[2]
        filename = sys.argv[3]
    except ValueError:
        exit_with_error()

    data_matrix = read_data_file(filename)

    if(len(data_matrix) <= k or len(data_matrix) == 0 or (k <= 1 and goal not in ('norm', 'ddg', 'sym'))):
        exit_with_error()

    if goal == "sym":
        sym(data_matrix)
    elif goal == "ddg":
        ddg(data_matrix)
    elif goal == "norm":
        norm(data_matrix)
    elif goal == "symnmf":
        display_matrix(symnmf_algo(data_matrix, k, BETA, EPSILON, MAX_ITERS))
        
    else:
        exit_with_error()

if __name__ == "__main__":
    main()