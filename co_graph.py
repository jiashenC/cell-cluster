import numpy as np
import distance as dis
from collections import deque
from scipy.stats import hypergeom
import igraph as ig
import louvain


class CoGraph:
    # data should be in numpy.ndarray format
    def __init__(self, data):
        self.pq = []
        self.data = data
        self.matrix = {}
        self.graph = None

    def co_test(self, i, j):
        left = self.data[i].astype(bool)
        right = self.data[j].astype(bool)
        k = np.count_nonzero(np.bitwise_and(left, right))
        prb = hypergeom.cdf(k, len(left), np.count_nonzero(left), np.count_nonzero(right))
        return 1 - prb

    def build_graph(self, threshold=0.05, jaccard=False, jaccard_threshold=0.5):
        self.graph = Graph()
        for i in range(0, len(self.data)):
            for j in range(i + 1, len(self.data)):
                co_score = self.co_test(i, j)
                if co_score > threshold:
                    self.matrix.update({(i, j): co_score})
                    self.graph.add(i, j)
        if jaccard:
            self.jaccard_preprocess(jaccard_threshold)

    # use BFS to give jaccard score to each pair and non-direct connected pairs
    def jaccard_preprocess(self, threshold):
        # track which two pairs have been tested
        test_list = []
        queue = deque([])
        queue.append(self.graph.vertices[0])
        while len(queue) != 0:
            cur = queue.popleft()
            neighbors = self.graph.vertices_matrix[cur]
            for node in neighbors:
                if (node, cur) in test_list or (cur, node) in test_list:
                    continue

                score = dis.jaccard(self.graph.vertices_matrix[cur], self.graph.vertices_matrix[node])
                if score < threshold:
                    self.graph.vertices_matrix[cur].remove(node)
                    self.graph.vertices_matrix[node].remove(cur)
                    if (cur, node) in self.graph.edges:
                        self.graph.edges.remove((cur, node))
                    else:
                        self.graph.edges.remove((node, cur))

                second_neighbors = self.graph.vertices_matrix[node]
                for second_node in second_neighbors:
                    if (cur, second_node) in test_list or (second_node, cur) in test_list or (
                            cur, second_node) in self.graph.edges or (second_node, cur) in self.graph.edges:
                        continue
                    test_list.append((cur, second_node))
                    second_score = dis.jaccard(self.graph.vertices_matrix[cur], self.graph.vertices_matrix[second_node])
                    if second_score > threshold:
                        self.graph.vertices_matrix[cur].append(second_node)
                        self.graph.vertices_matrix[second_node].append(cur)
                        self.graph.edges.append((cur, second_node))

        for vertex in self.graph.vertices:
            if len(self.graph.vertices_matrix[vertex]) < 1:
                self.graph.vertices.remove(vertex)

    def find_partition(self):
        g = ig.Graph(self.graph.edges)
        print g
        weights = []
        for u, v in self.graph.edges:
            if (u, v) in self.matrix:
                weights.append(self.matrix[(u, v)])
            elif (v, u) in self.matrix:
                weights.append(self.matrix[(v, u)])
            else:
                weights.append(self.co_test(u, v))
        g.es['weight'] = weights
        parts = louvain.find_partition(g, method='Modularity', weight='weight')
        print parts


class Graph:
    def __init__(self):
        # [y, x, ...]
        self.vertices = []
        # [(u,v), (x, y) ... ]
        self.edges = []
        # u -> [x, y .. ]
        self.vertices_matrix = {}

    def add(self, i, j):
        if (i, j) in self.edges or (j, i) in self.edges:
            return

        if i not in self.vertices_matrix:
            self.vertices.append(i)
            self.vertices_matrix.update({i: []})
        self.vertices_matrix[i].append(j)

        if j not in self.vertices_matrix:
            self.vertices.append(j)
            self.vertices_matrix.update({j: []})
        self.vertices_matrix[j].append(i)

        self.edges.append((i, j))
