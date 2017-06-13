import h5py
import numpy as np


class DataReader:
    def __init__(self, path='../data/raw_retina_data.mat', debug=False, debug_data=None):
        if not debug:
            self.file = h5py.File(path)
            self.original_data = self.file['data']
            self.data = None
        else:
            self.original_data = debug_data
            self.data = None

    """
        :param original specifies return original data or not
        :return data in a format that each entry is a cell
    """

    def get_cell_data(self, original=False):
        return self.original_data if original else self.data

    """
    :param original specifies return original data or not
    :return data in a format that each entry is a gene
    """

    def get_gene_data(self, original=False):
        return np.array(self.original_data).transpose() if original else np.array(self.data).transpose()

    def pre_process_cell(self, threshold=900, original=False):
        data = self.get_cell_data(original)
        tb = []
        for item in data:
            count = np.count_nonzero(item)
            if count > threshold:
                tb.append(item)
        self.data = np.array(tb)

    def pre_process_gene(self, threshold=900, original=False):
        data = self.get_gene_data(original)
        tb = []
        for item in data:
            count = np.count_nonzero(item)
            if count > threshold:
                tb.append(item)
        self.data = np.array(tb).transpose()
