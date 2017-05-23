import heapq as h
import scipy.stats as ss
import numpy as np
from sklearn.metrics import jaccard_similarity_score


class ChisqGraph:
    # data should be in numpy.ndarray format
    def __init__(self, data):
        self.pq = []
        self.data = data
        tb = []
        for i in range(0, len(data)):
            for j in range(i + 1, len(data)):
                if (i, j) not in tb or (j, i) not in tb:
                    h.heappush(self.pq, (-1.0 * self.chi_square_prob(i, j), i, j))
                    tb.append((i, j))

    def chi_square_prob(self, i, j):
        obs = np.array(self.data[i], self.data[j])
        chisq, chisq_prob = ss.chisquare(obs, axis=None)
        return chisq_prob

    def build_graph(self, threshold=0.05, jaccard=False, jaccard_threshold=0):
        graph = Graph()
        while True:
            chisq, u, v = h.heappop(self.pq)
            if abs(chisq) < threshold:
                break
            graph.add(u, v, chisq)
        if jaccard:
            self.jaccard_preprocess(graph, jaccard_threshold)

    def jaccard_preprocess(self, graph, threshold):
        for i, dist, j in graph.edges:
            jc = jaccard_similarity_score(self.data[i], self.data[j])
            if jc < threshold:
                graph.edges.remove((i, dist, j))
                graph.matrix[i].remove((j, dist))
                graph.matrix.pop(i, None)
                graph.matrix[j].remove((i, dist))
                graph.matrix.pop(j, None)
        # TODO: use further jaccard to connect not directly connected components

class Graph:
    def __init__(self):
        # vertices list
        self.vertices = []
        # edges list
        self.edges = []
        # graph matrix v -> (u, dist)
        self.matrix = {}

    def add(self, u, v, dist):
        if (u, dist, v) in self.edges:
            return
        else:
            if v not in self.vertices:
                self.vertices.append(v)
                self.matrix.update({v: []})
            self.matrix[v].append((u, dist))
            if u not in self.vertices:
                self.vertices.append(u)
                self.matrix.update({u: []})
            self.matrix[u].append((v, dist))
