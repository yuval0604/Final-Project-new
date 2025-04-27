import math
import sys 

def euclidean_distance(vec1, vec2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(vec1, vec2)))

def assign_clusters(data_points, centroids):
    labels = []
    for point in data_points:
        min_idx = 0
        min_dist = euclidean_distance(point, centroids[0])
        for i in range(1, len(centroids)):
            dist = euclidean_distance(point, centroids[i])
            if dist < min_dist:
                min_dist = dist
                min_idx = i
        labels.append(min_idx)
    return labels

def update_centroids(data_points, labels, k):
    d = len(data_points[0])
    centroids = [[0.0] * d for i in range(k)]
    cnt = [0] * k

    for point, label in zip(data_points, labels):
        for i in range(d):
            centroids[label][i] += point[i]
        cnt[label] += 1

    for i in range(k):
        if cnt[i] > 0:
            centroids[i] = [x / cnt[i] for x in centroids[i]]
    return centroids

def kmeans(data_points, k, max_iter=300, epsilon=1e-4):
    centroids = [point[:] for point in data_points[:k]]

    for i in range(max_iter):
        labels = assign_clusters(data_points, centroids)
        new_centroids = update_centroids(data_points, labels, k)

        converged = all(euclidean_distance(centroids[i], new_centroids[i]) < epsilon for i in range(k))

        centroids = new_centroids
        if converged:
            break

    return centroids, labels

def read_data(filename):
    try:
        with open(filename, 'r') as f:
            data = [list(map(float, line.strip().split(','))) for line in f if line.strip()]
        return data
    except:
        print("An Error Has Occurred", file=sys.stderr)
        sys.exit(1)

def parse_args():
    if not (3 <= len(sys.argv) <= 4):
        print("An Error Has Occurred")
        sys.exit(1)

    try:
        k = int(sys.argv[1])
    except:
        print("Invalid number of clusters!")
        sys.exit(1)
    
    try:
        max_iter = int(sys.argv[2]) if len(sys.argv) == 4 else 300
    except:
        print("Invalid maximum iteration!")
        sys.exit(1)

    input_file = sys.argv[3] if len(sys.argv) == 4 else sys.argv[2]

    if max_iter <= 1 or max_iter >= 1000:
        print("Invalid maximum iteration!")
        sys.exit(1)
    
    if k <= 1:
        print("Invalid number of clusters!")
        sys.exit(1)

    return k, max_iter, input_file

if __name__ == "__main__":
    k, max_iter, input_file = parse_args()
    data_points = read_data(input_file)

    if k >= len(data_points):
        print("Invalid number of clusters!")
        sys.exit(1)

    centroids, labels = kmeans(data_points, k, max_iter)

    for c in centroids:
        print(','.join(f'{coord:.4f}' for coord in c))