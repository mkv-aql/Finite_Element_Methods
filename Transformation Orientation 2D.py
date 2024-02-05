__author__ = 'mkv-aql'
import numpy as np
from fractions import Fraction
import sympy as sp

#Note: * is not multiplication in matrix, use x.dot(y)

#Decimal point
np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)}, precision = 2)

# Angle of the element must be in radians
alpha = 150/360*2*np.pi

K_Matrix = np.array([[1, -1],
                    [-1, 1]])

Transformation = np.array([[np.cos(alpha), np.sin(alpha), 0, 0],
                          [0, 0, np.cos(alpha), np.sin(alpha)]])

print(Transformation)
# print(Transformation.transpose())

Transformation_transpose = Transformation.transpose()

#K = T tranpose * Stiffness matrix * T
Stiffness_Matrix_Transformed = Transformation_transpose.dot(K_Matrix.dot(Transformation))

print("\n", Stiffness_Matrix_Transformed)

