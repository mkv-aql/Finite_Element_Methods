__author__ = 'mkv-aql'

import numpy as np

# Parameters
E = 200e9
v = 0.3


# Formula
c_matrix = np.array([[1, v, 0],
                     [v, 1, 0],
                     [0, 0, (1-v)/2]]) * E / (1 - v**2)

# round off to 2 decimal places
c_matrix = np.round(c_matrix, 2)

print(c_matrix)

def Material_Matrix(E, v):
    c_matrix = np.array([[1, v, 0],
                     [v, 1, 0],
                     [0, 0, (1-v)/2]]) * E / (1 - v**2)
    c_matrix = np.round(c_matrix, 2)
    return c_matrix