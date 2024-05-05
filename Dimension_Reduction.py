__author__ = 'mkv-aql'

import numpy as np
from scipy.interpolate import interp2d

# Create a 6x6 matrix
matrix_6x6 = np.array([
    [600, -600, 0, 0, 0, 0],
    [-600, 1200, -600, 0, 0, 0],
    [0, -600, 1200, -600, 0, 0],
    [0, 0, -600, 1000, -400, 0],
    [0, 0, 0, -400, 800, -400],
    [0, 0, 0, 0, -400, 400]
])

# Create interpolation function
x = y = np.arange(6)
interp_func = interp2d(x, y, matrix_6x6, kind='linear')

# New coordinates for the 5x5 matrix
new_x = new_y = np.linspace(0, 5, 5)

# Use the interpolation function to create the 5x5 matrix
matrix_5x5 = interp_func(new_x, new_y)

print(matrix_5x5)
