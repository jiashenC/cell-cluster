from data_reader import DataReader
from co_graph import CoGraph
import time
import numpy as np
from scipy.stats import hypergeom, chisquare
import validation
import log as lg
import igraph as ig
import louvain

# arr = np.array([[1, 2, 3, 4], [2, 3, 4, 5]])
# print arr.shape
# print arr[0]
# arr = arr.transpose()
# print arr.shape
# print arr[0]

# arr = np.array([1., 2., 3., 4., 0., 0.])
# arr = arr.astype(bool)
# arr1 = np.array([0., 0., 3., 4., 1., 2.])
# arr1 = arr1.astype(bool)
# print np.bitwise_and(arr, arr1)

log = lg.Log()
data = DataReader()
data.pre_process_cell(original=True)
data.pre_process_gene(threshold=6000, original=False)
graph = CoGraph(data.get_gene_data())
graph.build_graph(threshold=0.5, jaccard=True, jaccard_threshold=0.5)
validation.valid_graph(graph.graph.vertices, graph.graph.edges, graph.graph.vertices_matrix)
log.update()
graph.find_partition()
log.time('community detection algorithm')
