import h5py
import numpy as np


class DataReader:
    def __init__(self, path='../data/raw_retina_data.mat', debug=False, debug_data=None):
        if not debug:
            self.file = h5py.File(path)
            self.data = None
            group_keys = self.file.keys()
            for k in range(0, len(group_keys)):
                self.data = np.array(self.file[group_keys[k]]) if str(
                    group_keys[k]) == 'data' or self.data is not None else None
        else:
            self.data = debug_data

    # return row data in one array
    def get_row_data(self):
        return self.data

    # return column data in one array
    def get_col_data(self):
        return self.data.transpose()


