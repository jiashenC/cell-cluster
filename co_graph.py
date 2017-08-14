import numpy as np
import distance as dis
from collections import deque
from scipy.stats import hypergeom
from scipy.spatial.distance import euclidean
import igraph as ig
import louvain
from log import output
import Queue
import threading


class CoGraph:
    # data should be in numpy.ndarray format
    def __init__(self, data):
        self.pq = []
        self.data = data
        self.data_length = len(data)
        self.matrix = {}
        self.graph = None
        self.parts = None
        self.queue = Queue.Queue()
        self.rqueue = Queue.Queue()
        self.test_threshold = 0
        self.finish = False

    def co_test_single(self, i, j):
        left = self.data[i].astype(bool)
        right = self.data[j].astype(bool)
        k = np.count_nonzero(np.bitwise_and(left, right))
        prb = hypergeom.cdf(k, len(left), np.count_nonzero(left), np.count_nonzero(right))
        return 1 - prb

    def eu_test_single(self, i, j):
        return euclidean(self.data[i], self.data[j])

    def co_test(self):
        while True:
            i, j = self.queue.get()
            left = self.data[i].astype(bool)
            right = self.data[j].astype(bool)
            k = np.count_nonzero(np.bitwise_and(left, right))
            prb = hypergeom.cdf(k, len(left), np.count_nonzero(left), np.count_nonzero(right))
            if 1 - prb < self.test_threshold:
                self.rqueue.put((i, j))
            if i == self.data_length - 2:
                self.finish = True
                self.rqueue.put((-1, -1))
                break

    def eu_test(self):
        while True:
            i, j = self.queue.get()
            if euclidean(self.data[i], self.data[j]) < self.test_threshold:
                self.rqueue.put((i, j))
            if i == self.data_length - 2:
                self.finish = True
                self.rqueue.put((-1, -1))
                break

    def build_graph(self, threshold=0.05, jaccard=False, jaccard_threshold=0.5, mode='hypergeometry'):
        self.test_threshold = threshold
        if mode == 'hypergeometry':
            threading.Thread(target=self.co_test).start()
        else:
            threading.Thread(target=self.eu_test).start()

        threading.Thread(target=self.add_test).start()

        self.graph = Graph()
        while not self.finish:
            i, j = self.rqueue.get()
            if i < 0:
                break
            self.graph.add(i, j)

        output('after ' + mode + ' edge: ', len(self.graph.edges))
        if jaccard:
            self.jaccard_preprocess(jaccard_threshold)

    def add_test(self):
        for i in range(0, self.data_length):
            for j in range(i + 1, self.data_length):
                self.queue.put((i, j))

    # use BFS to give jaccard score to each pair and non-direct connected pairs
    def jaccard_preprocess(self, threshold):
        # track which two pairs have been tested
        test_edge_list = set([])
        test_vertices_list = set([])
        queue = deque([])
        count = 0
        for x in self.graph.vertices:
            queue.append(x)
            break
        while len(queue) != 0:
            cur = queue.popleft()
            if cur in test_vertices_list:
                continue

            count += 1
            test_vertices_list.add(cur)
            neighbors = self.graph.vertices_matrix[cur]
            for node in neighbors:
                if (node, cur) in test_edge_list or (cur, node) in test_edge_list:
                    continue

                test_edge_list.add((cur, node))
                if node not in test_vertices_list:
                    queue.append(node)
                score = dis.jaccard(self.graph.vertices_matrix[cur], self.graph.vertices_matrix[node])
                self.matrix.update({(cur, node): score})
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
                    self.matrix.update({(cur, second_node): second_score})
                    if second_score > threshold:
                        self.graph.vertices_matrix[cur].append(second_node)
                        self.graph.vertices_matrix[second_node].append(cur)
                        self.graph.edges.add((cur, second_node))

        for vertex in self.graph.vertices:
            if len(self.graph.vertices_matrix[vertex]) < 1:
                self.graph.vertice.remove(vertex)

        output('after jaccard', len(self.graph.edges))

    def find_partition(self, weight=True, mode='hypergeometry'):
        g = ig.Graph(list(self.graph.edges))
        # use hyper geometry test as edge weights
        weights = []
        for u, v in self.graph.edges:
            if (u, v) in self.matrix:
                weights.append(self.matrix[(u, v)])
            elif (v, u) in self.matrix:
                weights.append(self.matrix[(v, u)])
            else:
                score = self.co_test_single(u, v) if mode == 'hypergeometry' else self.eu_test_single(u, v)
                weights.append(score)
        if weight:
            g.es['weight'] = weights
            self.parts = louvain.find_partition(g, method='Modularity', weight='weight')
        else:
            self.parts = louvain.find_partition(g, method='Modularity')


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
