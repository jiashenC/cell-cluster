from scipy.stats import hypergeom
import numpy as np


def valid_graph(vertices, edges, vertices_matrix):
    for v in vertices:
        if v not in vertices_matrix:
            raise Exception('vertices list and vertices_matrix is not equal')

    for u, v in edges:
        if u not in vertices_matrix or v not in vertices_matrix[u]:
            raise Exception('vertices list and vertices_matrix is not equal')

    for u in vertices_matrix:
        for v in vertices_matrix[u]:
            if (u, v) not in edges and (v, u) not in edges:
                raise Exception('edge list and vertices_matrix is not equal')


def valid_covalue(data, parts):
    count = 0
    total = 0
    for i in range(0, len(parts) - 1):
        for j in range(i + 1, len(parts)):
            for partone in range(0, len(parts[i])):
                for parttwo in range(0, len(parts[j])):
                    left_i = parts[i][partone]
                    right_i = parts[j][parttwo]
                    left = data[left_i].astype(bool)
                    right = data[right_i].astype(bool)
                    k = np.count_nonzero(np.bitwise_and(left, right))
                    prb = hypergeom.cdf(k, len(left), np.count_nonzero(left), np.count_nonzero(right))
                    if 1 - prb < 0.05:
                        count += 1
                    total += 1
    print 'false positive %d' % ((1.0 * count) / total)

    count = 0
    total = 0
    for i in range(0, len(parts)):
        for k in range(0, len(parts[i])):
            for j in range(k, len(parts[i])):
                left_i = parts[i][k]
                right_i = parts[i][j]
                left = data[left_i].astype(bool)
                right = data[right_i].astype(bool)
                k = np.count_nonzero(np.bitwise_and(left, right))
                prb = hypergeom.cdf(k, len(left), np.count_nonzero(left), np.count_nonzero(right))
                if 1 - prb > 0.05:
                    count += 1
                total += 1
    print 'true positive %d' % (1.0 * count / total)
