__author__ = 'mkv-aql'

import numpy as np

import Function_Module as sm

print(sm.Stiffness_Matrix_1(180, 180))
print(sm.Stiffness_Matrix_1(0, 0))
print(sm.Stiffness_Matrix_2(180))
print(sm.Stiffness_Matrix_2(0))


a = np.array([[1200, -600, 0, 0, 0], [-600, 1200, -600, 0, 0], [0, -600, 1000, -400, 0], [0, 0, -400, 800, -400], [0, 0, 0, -400, 400]])
x = 25/3
b = np.array([x+x, x+x, x, 0, 25])

print(np.linalg.solve(a, b))
