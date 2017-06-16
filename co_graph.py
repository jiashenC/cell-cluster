import numpy as np
import distance as dis
from collections import deque
from scipy.stats import hypergeom
import igraph as ig
import louvain
import log


class CoGraph:
    # data should be in numpy.ndarray format
    def __init__(self, data):
        self.pq = []
        self.data = data
        # self.matrix = {}
        self.graph = None

    def co_test(self, i, j):
        left = self.data[i].astype(bool)
        right = self.data[j].astype(bool)
        k = np.count_nonzero(np.bitwise_and(left, right))
        prb = hypergeom.cdf(k, len(left), np.count_nonzero(left), np.count_nonzero(right))
        return 1 - prb

    def build_graph(self, threshold=0.05, jaccard=False, jaccard_threshold=0.5):
        self.graph = Graph()
        lg = log.Log()
        for i in range(0, len(self.data)):
            for j in range(i + 1, len(self.data)):
                co_score = self.co_test(i, j)
                if co_score < threshold:
                    # self.matrix.update({(i, j): co_score})
                    self.graph.add(i, j)
        lg.time('finish hyper geometry')
        if jaccard:
            self.jaccard_preprocess(jaccard_threshold)
        lg.time('finish jaccard')

    # use BFS to give jaccard score to each pair and non-direct connected pairs
    def jaccard_preprocess(self, threshold):
        # track which two pairs have been tested
        test_edge_list = set([])
        test_vertices_list = set([])
        queue = deque([])
        for x in self.graph.vertices:
            queue.append(x)
            break
        while len(queue) != 0:
            cur = queue.popleft()
            if cur in test_vertices_list:
                continue

            test_vertices_list.add(cur)
            neighbors = self.graph.vertices_matrix[cur]
            for node in neighbors:
                if (node, cur) in test_edge_list or (cur, node) in test_edge_list:
                    continue

                test_edge_list.add((cur, node))
                queue.append(node)
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
                    if (cur, second_node) in test_edge_list or (second_node, cur) in test_edge_list or (
                            cur, second_node) in self.graph.edges or (second_node, cur) in self.graph.edges:
                        continue
                    test_edge_list.add((cur, second_node))
                    second_score = dis.jaccard(self.graph.vertices_matrix[cur], self.graph.vertices_matrix[second_node])
                    if second_score > threshold:
                        self.graph.vertices_matrix[cur].append(second_node)
                        self.graph.vertices_matrix[second_node].append(cur)
                        self.graph.edges.add((cur, second_node))

        for vertex in self.graph.vertices:
            if len(self.graph.vertices_matrix[vertex]) < 1:
                self.graph.vertice.remove(vertex)

    def find_partition(self):
        g = ig.Graph(list(self.graph.edges))
        # use hyper geometry test as edge weights
        # weights = []
        # for u, v in self.graph.edges:
        #     if (u, v) in self.matrix:
        #         weights.append(self.matrix[(u, v)])
        #     elif (v, u) in self.matrix:
        #         weights.append(self.matrix[(v, u)])
        #     else:
        #         weights.append(self.co_test(u, v))
        # g.es['weight'] = weights
        parts = louvain.find_partition(g, method='Modularity')
        print parts


class Graph:
    def __init__(self):
        # [y, x, ...]
        self.vertices = set([])
        # [(u,v), (x, y) ... ]
        self.edges = set([])
        # u -> [x, y .. ]
        self.vertices_matrix = {}

    def add(self, i, j):
        if (i, j) in self.edges or (j, i) in self.edges:
            return

        if i not in self.vertices_matrix:
            self.vertices.add(i)
            self.vertices_matrix.update({i: []})
        self.vertices_matrix[i].append(j)

        if j not in self.vertices_matrix:
            self.vertices.add(j)
            self.vertices_matrix.update({j: []})
        self.vertices_matrix[j].append(i)

        self.edges.add((i, j))
