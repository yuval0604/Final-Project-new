import sys
import numpy as np
from sklearn.metrics import silhouette_score
import kmeans
import symnmf


ERROR_MSG = "An Error Has Occurred"
BETA = 0.5
EPSILON = 1e-4
MAX_ITERS = 300

def read_data_file(filepath):
    """Reads a text file into a representing a matrix."""
    return symnmf.read_data_file(filepath)

def kmeans_analysis(data, k, max_iter, epsilon):
    """
    Runs KMeans clustering and returns the labels.

    :param data: input data as numpy array
    :param k: number of clusters
    :param max_iter: max number of iterations
    :param epsilon: convergence threshold
    :return: labels (list of cluster assignments)
    """
    centroids, labels = kmeans.kmeans(data.tolist(), k, max_iter=max_iter, epsilon=epsilon)
    return labels

def symnmf_analysis(data_matrix, k_value, beta, epsilon, max_iter):
    data_matrix = [[float(x) for x in row] for row in data_matrix.tolist()]  
    H_matrix = symnmf.symnmf_algo(data_matrix, k_value, beta, epsilon, max_iter)
    return np.argmax(H_matrix, axis=1)
    
def main():
    if len(sys.argv) != 3:
        print(ERROR_MSG)
        sys.exit(1)

    try:
        path = sys.argv[2]
        k = int(sys.argv[1])
    except Exception:
        print(ERROR_MSG)
        sys.exit(1)

    data = read_data_file(path)
    data_np = np.array(data)
    try :
        kmeans_labels = kmeans_analysis(data_np, k, MAX_ITERS, EPSILON)
        kmeans_score = silhouette_score(data_np, kmeans_labels)
        symnmf_lables = symnmf_analysis(data_np, k, BETA, EPSILON, MAX_ITERS)
        symnmf_score = silhouette_score(data_np, symnmf_lables)
        print(f"nmf: {symnmf_score:.4f}")
        print(f"kmeans: {kmeans_score:.4f}")
    except:
        print(ERROR_MSG)
        sys.exit(1)
    

if __name__ == "__main__":
    main()