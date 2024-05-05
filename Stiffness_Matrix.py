__author__ = 'mkv-aql'

import numpy as np
from fractions import Fraction
import sympy as sp

# Angle of the element must be in radians
alpha = 45/360*2*np.pi
beta = 0/360*2*np.pi

# Stiffness matrix for a single element
stiffness_matrix = np.array([[np.cos(alpha)**2, np.sin(alpha)*np.cos(alpha), -np.cos(alpha)*np.cos(beta), -np.sin(beta)*np.cos(alpha)],
                                [np.sin(alpha)*np.cos(alpha), np.sin(alpha)**2, -np.sin(alpha)*np.cos(beta), -np.sin(alpha)*np.sin(beta)],
                                [-np.cos(alpha)*np.cos(beta), -np.sin(alpha)*np.cos(beta), np.cos(beta)**2, np.sin(beta)*np.cos(beta)],
                                [-np.sin(beta)*np.cos(alpha), -np.sin(alpha)*np.sin(beta), np.sin(beta)*np.cos(beta), np.sin(beta)**2]])

# print(stiffness_matrix)
# print(stiffness_matrix*4)

'''
dec_value = 0.666666666
rational_value = Fraction(dec_value).limit_denominator()
print(rational_value)

# loop through the matrix and convert the decimal values to rational values
for i in range(0,4):
    for j in range(0,4):
        rational_value = Fraction(stiffness_matrix[i][j]).limit_denominator(10000)
        #print(rational_value)
        stiffness_matrix[i][j] = rational_value

# print the stiffness matrix
fraction_stiffness_matrix = np.vectorize(Fraction)(stiffness_matrix)
print(fraction_stiffness_matrix)
'''
'''
from fractions import Fraction
import math


def decimal_to_fraction_with_sqrt(decimal):
    if decimal < 0:
        return f"-√{abs(decimal):.4f}"

    sqrt_value = math.isqrt(int(decimal * 1e8))
    fraction = Fraction(sqrt_value, 1e4)

    return f"√{decimal:.4f} = {fraction}"


# Example decimal value
decimal_value = 0.7071

fraction_representation = decimal_to_fraction_with_sqrt(decimal_value)
print(fraction_representation)
'''

def Stiffness_Matrix_1(alpha, beta):
    alpha = alpha / 360 * 2 * np.pi
    beta = beta / 360 * 2 * np.pi

    stiffness_matrix = np.array([[np.cos(alpha) ** 2, np.sin(alpha) * np.cos(alpha), -np.cos(alpha) * np.cos(beta),
                                  -np.sin(beta) * np.cos(alpha)],
                                 [np.sin(alpha) * np.cos(alpha), np.sin(alpha) ** 2, -np.sin(alpha) * np.cos(beta),
                                  -np.sin(alpha) * np.sin(beta)],
                                 [-np.cos(alpha) * np.cos(beta), -np.sin(alpha) * np.cos(beta), np.cos(beta) ** 2,
                                  np.sin(beta) * np.cos(beta)],
                                 [-np.sin(beta) * np.cos(alpha), -np.sin(alpha) * np.sin(beta),
                                  np.sin(beta) * np.cos(beta), np.sin(beta) ** 2]])
    return stiffness_matrix

def Stiffness_Matrix_Mapping_2D_Single(b_matrix, c_matrix, det_jacobian, weight):
    # B transpose
    b_matrix_transpose = b_matrix.transpose()

    # Stiffness Matrix Mapping
    stiffness = b_matrix_transpose * c_matrix * b_matrix * det_jacobian * weight

    return stiffness

