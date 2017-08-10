import h5py
import numpy as np


class DataReader:
    def __init__(self, path='../data/raw_retina_data.mat', debug=False, debug_data=None):
        if not debug:
            self.file = h5py.File(path)
            self.original_data = self.file['data']
            self.data = None
            self.cell_lookup = {}
            self.gene_lookup = {}
        else:
            self.original_data = debug_data
            self.data = None

    def get_cell_data(self, original=False):
        """
            :param original specifies return original data or not
            :return data in a format that each entry is a cell
        """
        return self.original_data if original else self.data

    def get_gene_data(self, original=False):
        """
            :param original specifies return original data or not
            :return data in a format that each entry is a gene
        """
        return np.array(self.original_data).transpose() if original else np.array(self.data).transpose()

    def pre_process_cell(self, threshold=900, original=False):
        data = self.get_cell_data(original)
        n_index = 0
        pre_index = 0
        tb = []
        for item in data:
            count = np.count_nonzero(item)
            if (threshold > 1 and count > threshold) or (threshold < 1 and count * 1.0 / len(item) > 1 - threshold):
                tb.append(item)
                self.cell_lookup.update({n_index: pre_index})
                n_index += 1
                pre_index += 1
            else:
                pre_index += 1
        self.data = np.array(tb)

    def pre_process_gene(self, threshold=900, original=False):
        data = self.get_gene_data(original)
        n_index = 0
        pre_index = 0
        tb = []
        for item in data:
            count = np.count_nonzero(item)
            if (threshold > 1 and count > threshold) or (threshold < 1 and count * 1.0 / len(item) > 1 - threshold):
                tb.append(item)
                self.gene_lookup.update({n_index: pre_index})
                n_index += 1
                pre_index += 1
            else:
                pre_index += 1
        self.data = np.array(tb).transpose()
