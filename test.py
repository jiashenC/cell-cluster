from data_reader import DataReader
import numpy as np

arr = np.array([[1, 2, 3], [1, 2, 3], [1, 2, 3], [1, 2, 3]])
data = DataReader('hello', True, arr)
print data.get_col_data()
print data.get_row_data()
