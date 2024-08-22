__author__ = 'mkv-aql'
import math

from sympy import Symbol, cos, factor, symbols, integrate

import numpy as np

#weights
wg = 1

# Coordinates of the nodes, bottom left, bottom right, top right, top left
# cor = [[1, 1], [6, 1], [6, 4], [1, 4]] # Parameter
cor = [[0, 0], [1.5, 0], [1.5, 2/3], [0, 2/3]] # Parameter

# Area of the element
area = (cor[1][0] - cor[0][0]) * (cor[2][1] - cor[0][1]) # Parameter
print(f'\nArea of the element: {area}')

# stiffness k
stiffness_k = 1 # Parameter
# stiffness_k = 500.0556731 # Parameter

# Linear shape functions for quadrilateral element, h = high, l = low
x, y, x_l, x_h, y_l, y_h, inv_area, k = symbols('x, y, x_l, x_h, y_l, y_h, 1/A, k')
h1 = (inv_area)*(x_h - x)*(y_h - y)
h2 = (inv_area)*(x - x_l)*(y_h - y)
h3 = (inv_area)*(x - x_l)*(y - y_l)
h4 = (inv_area)*(x_h - x)*(y - y_l)

print("Linear shape functions formulas: ")
print("h1 = ", h1)
print("h2 = ", h2)
print("h3 = ", h3)
print("h4 = ", h4)

# Bilienar shape functions for quadrilateral element, h = high, l = low
dh1_dx = h1.diff(x)
dh1_dy = h1.diff(y)
dh2_dx = h2.diff(x)
dh2_dy = h2.diff(y)
dh3_dx = h3.diff(x)
dh3_dy = h3.diff(y)
dh4_dx = h4.diff(x)
dh4_dy = h4.diff(y)

print("\nBilinear shape functions formulas: ")
print("dh1_dx = ", dh1_dx)
print("dh1_dy = ", dh1_dy)
print("dh2_dx = ", dh2_dx)
print("dh2_dy = ", dh2_dy)
print("dh3_dx = ", dh3_dx)
print("dh3_dy = ", dh3_dy)
print("dh4_dx = ", dh4_dx)
print("dh4_dy = ", dh4_dy)
print(f'\nNOTE: the formulas seem to be wrong, but it is actually correct, - and + are simplified')

# B Matrix
b1 = np.array([dh1_dx, dh1_dy])
b2 = np.array([dh2_dx, dh2_dy])
b3 = np.array([dh3_dx, dh3_dy])
b4 = np.array([dh4_dx, dh4_dy])

B = np.array([b1[0], b2[0], b3[0], b4[0],
              b1[1], b2[1], b3[1], b4[1]])
print(f"\nB Matrix: \n{B}")


# b2 x b3
result = np.dot(b2, b3)
print(f"\nb2 x b3 = {result}")
expanded_result = result.expand()
factor_result = result.factor()
print(f'\nExpanded: {expanded_result}')
print(f'\nFactored: {factor_result}')

# Sub stiffness matrix, k = stiffness, s = sub stiffness, inverse of area already calculated into each b
s11 = k * np.dot(b1, b1)
s12 = k * np.dot(b1, b2)
s13 = k * np.dot(b1, b3)
s14 = k * np.dot(b1, b4)
s21 = k * np.dot(b2, b1)
s22 = k * np.dot(b2, b2)
s23 = k * np.dot(b2, b3)
s24 = k * np.dot(b2, b4)
s31 = k * np.dot(b3, b1)
s32 = k * np.dot(b3, b2)
s33 = k * np.dot(b3, b3)
s34 = k * np.dot(b3, b4)
s41 = k * np.dot(b4, b1)
s42 = k * np.dot(b4, b2)
s43 = k * np.dot(b4, b3)
s44 = k * np.dot(b4, b4)

print(f'\nK23 formula: {s23}')
print(f'K11 formula: {s11.factor()}')

# insert values
# s23 = s23.subs({inv_area: 1/area, x_l: cor[0][0], x_h: cor[1][0], y_l: cor[0][1], y_h: cor[2][1]})

# intergrals of the sub stiffness matrix
# s23 = integrate(s23, (x, cor[0][0], cor[1][0]), (y, cor[0][1], cor[2][1]))

print(f'\nSub stiffness matrix: K23 \n{s23}')

K = [s11, s21, s31, s41,
    s12, s22, s32, s42,
    s13, s23, s33, s43,
    s14, s24, s34, s44]

# insert values for all sub stiffness matrix, k=1 for simplicity
K = [i.subs({inv_area: 1/area, x_l: cor[0][0], x_h: cor[1][0], y_l: cor[0][1], y_h: cor[2][1], k: stiffness_k}) for i in K]

# intergrals of the stiffness matrix
K = [integrate(i, (x, cor[0][0], cor[1][0]), (y, cor[0][1], cor[2][1])) for i in K]

# convert to 4 decimal places
K = [round(i, 4) for i in K]

K = np.array(K).reshape(4, 4)
print(f'\nStiffness matrix of coordinate: \n{cor} \nK: \n{K}')



