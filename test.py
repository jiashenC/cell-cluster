from data_reader import DataReader
from co_graph import CoGraph
import time
import numpy as np
import validation
import log as lg
import cluster_merge as cm
import igraph as ig
import louvain
from ZIFA import ZIFA
from ZIFA import block_ZIFA

# rd = DataReader()
# rd.pre_process_gene(threshold=0.8, original=True)
# data = np.log2(rd.get_cell_data(original=False) + 1)
# print data.shape
# Z, model_params = block_ZIFA.fitModel(data, 2)
# print Z

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

rd = DataReader()
rd.pre_process_cell(threshold=2700, original=True)
rd.pre_process_gene(threshold=1000, original=False)

g = CoGraph(rd.get_gene_data())
print g.data.shape
g.build_graph(threshold=0.001, jaccard=True, jaccard_threshold=0.5)

log.update()

g.find_partition()
# validation.valid_covalue(g.data, g.parts)
cm = cm.Cluster(g.data, g.parts)
cm.merge()

g2 = CoGraph(cm.parts)
print g2.data.shape
g2.build_graph(threshold=2, jaccard=False, mode='euclidean')
g2.find_partition()
print g2.parts

with open('../data/result.txt') as f:
    data = f.readlines()
    cluster = [(x.strip()).split()[1] for x in data]

re = np.zeros((39, len(g2.parts)))
col = 0
for part in g2.parts:
    for num in part:
        index = rd.cell_lookup[num]
        row = int(cluster[index])
        re[row - 1][col] += 1
    col += 1

print re
# row_sums = re.sum(axis=1)
# new_matrix = re / row_sums[:, np.newaxis]
# print new_matrix

log.time('community detection algorithm')
