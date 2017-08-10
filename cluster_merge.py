import numpy as np


class Cluster:
    def __init__(self, data, parts):
        self.data = data
        self.parts = parts

    def merge(self):
        tmp = []
        for k in range(0, len(self.parts)):
            array = []
            for index in range(0, len(self.parts[k])):
                array.append(self.data[index])
            arr = np.vstack(array)
            tmp.append(np.average(arr, axis=0, weights=arr.astype(bool)))
        self.parts = np.array(tmp).transpose()
