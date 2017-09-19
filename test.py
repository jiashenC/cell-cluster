from data_reader import DataReader
from co_graph import CoGraph
import numpy as np
import cluster_merge as cm
from log import output
import matplotlib.pyplot as plt
import scipy
from scipy.linalg import logm, expm

rd = DataReader()
rd.pre_process_cell(threshold=1850, original=True)
rd.pre_process_gene(threshold=2600, original=False)

genes = rd.get_gene_data()

collect = dict()

for x in range(len(genes)):
    for y in range(x + 1, len(genes)):
        X = genes[x]
        Y = genes[y]
        # matrix dimension
        N = 15
        matrix = [[0 for _ in range(N)] for _ in range(N)]
        for k in range(len(X)):
            xcord = X[k]
            ycord = Y[k]
            if xcord >= N or ycord >= N: continue
            matrix[int(xcord)][int(ycord)] += 1

        # for i in range(len(matrix[0])):
        #     matrix[0][i] = matrix[0][i] if matrix[0][i] != 0 else 1

        # from math import log
        # matrix = [[log(item if item != 0 else 1, 10) for item in row] for row in matrix]

        for k in range(len(matrix) - 1, -1, -1):
            result = ''
            for j in range(len(matrix[k])):
                result += '{:4}'.format(matrix[k][j])

        for k in range(1, N - 5):
            # for i in range(len(matrix[k])):
            #     matrix[k][i] = matrix[k][i] if matrix[k][i] != 0 else 1
            # print matrix[k]
            # print matrix[0]
            chi = False
            merge = [[0 for _ in range(N)] for _ in range(2)]
            for index in range(N):
                for i in range(k):
                    merge[0][index] += matrix[i][index]
                merge[0][index] = merge[0][index] / ((k + 1) * 1.0)
                if merge[0][index] == 0.0:
                    chi = True

            for index in range(N):
                for i in range(N - k):
                    merge[1][index] += matrix[i][index]
                merge[1][index] = merge[1][index] / ((20 - k) * 1.0)
                if merge[1][index] == 0.0:
                    chi = True

            if chi:
                continue

            chisq, p_value = scipy.stats.chisquare(merge[0], f_exp=merge[1])
            if p_value >= 0.05:
                collect.update({k: collect[k] + 1 if k in collect else 1})

print collect
