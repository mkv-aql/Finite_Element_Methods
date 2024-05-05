__author__ = 'mkv-aql'
import numpy as np
import pandas as pd

variable = ['1', '2', '3', '4']
matrix = np.array([[1, 2, 3, 4],
                   [5, 6, 7, 8],
                   [9, 10, 11, 12],
                   [13, 14, 15, 16]])

df = pd.DataFrame(matrix, columns=variable, index=variable)

print(df)

#reorder the matrix
df = df.reindex(columns=['4', '2', '3', '1'], index=['4', '2', '3', '1'])
print(df)