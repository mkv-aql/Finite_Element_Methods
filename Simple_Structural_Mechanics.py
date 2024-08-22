__author__ = 'mkv-aql'
import math

from sympy import Symbol, cos, factor, symbols, integrate

import numpy as np

# Parameters
# Coordinates of the nodes, bottom left, bottom right, top right, top left
cor = [[1, 1], [6, 1], [6, 4], [1, 4]] # Parameter
# cor = [[0, 0], [1.5, 0], [1.5, 2/3], [0, 2/3]] # Parameter

v_value = 0.3 # Parameter
E_value = 2.1e5 # Parameter

# Area of the element
area = (cor[1][0] - cor[0][0]) * (cor[2][1] - cor[0][1]) # Parameter
print(f'\nArea of the element: {area}')

# Linear shape functions for quadrilateral element, h = high, l = low
x, y, x_l, x_h, y_l, y_h, inv_area, v, E = symbols('x, y, x_l, x_h, y_l, y_h, 1/A, v, E')
h1 = (inv_area)*(x_h - x)*(y_h - y)
h2 = (inv_area)*(x - x_l)*(y_h - y)
h3 = (inv_area)*(x - x_l)*(y - y_l)
h4 = (inv_area)*(x_h - x)*(y - y_l)

C_1 = (E/(1-v**2))
C_2 = np.array([[1, v, 0],
                [v, 1, 0],
                [0, 0, (1-v)/2]])
print(f'\nC Matrix: \n{C_1} x \n{C_2}')

print("\nLinear shape functions formulas: ")
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
B1 = np.array([[dh1_dx, 0], [0, dh1_dy], [dh1_dy, dh1_dx]])
B2 = np.array([[dh2_dx, 0], [0, dh2_dy], [dh2_dy, dh2_dx]])
B3 = np.array([[dh3_dx, 0], [0, dh3_dy], [dh3_dy, dh3_dx]])
B4 = np.array([[dh4_dx, 0], [0, dh4_dy], [dh4_dy, dh4_dx]])
B1_Transpose = np.transpose(B1)
B2_Transpose = B2.transpose()
B3_Transpose = B3.transpose()
B4_Transpose = B4.transpose()

print(f"\nB1 Matrix: \n{B1}")
print(f"\nB2 Matrix: \n{B2}")
print(f"\nB3 Matrix: \n{B3}")
print(f"\nB4 Matrix: \n{B4}")

B = np.array([[B1[0][0], B1[0][1], B2[0][0], B2[0][1], B3[0][0], B3[0][1], B4[0][0], B4[0][1]],
              [B1[1][0], B1[1][1], B2[1][0], B2[1][1], B3[1][0], B3[1][1], B4[1][0], B4[1][1]],
              [B1[2][0], B1[2][1], B2[2][0], B2[2][1], B3[2][0], B3[2][1], B4[2][0], B4[2][1]]])
# print(f"\nB Matrix: \n{B}")

# B2_T x C x B1
result = np.dot(B2_Transpose, C_2)
result = np.dot(result, B3)
print(f"\nB2_T x C x B3 = \nSize: {result.shape}")
print(f'|{result[0][0]} | {result[0][1]} \n|{result[1][0]} | {result[1][1]}')
# expanded_result = result.expand()
# factor_result = result.factor()
# for i in result:
#     i[0].factor()
# print(f'\n|{result[0][0]} | {result[0][1]} \n|{result[1][0]} | {result[1][1]}')
# print(f'\nExpanded: {expanded_result}')
# print(f'\nFactored: {factor_result}')

# Sub stiffness matrix, k = stiffness, s = sub stiffness, inverse of area already calculated into each b
s11 = (np.dot(B1_Transpose, C_2)).dot(B1)
s12 = (np.dot(B2_Transpose, C_2)).dot(B1)
s13 = (np.dot(B2_Transpose, C_2)).dot(B1)
s14 = (np.dot(B4_Transpose, C_2)).dot(B1)
s21 = (np.dot(B1_Transpose, C_2)).dot(B2)
s22 = (np.dot(B2_Transpose, C_2)).dot(B2)
s23 = (np.dot(B3_Transpose, C_2)).dot(B2)
s24 = (np.dot(B4_Transpose, C_2)).dot(B2)
s31 = (np.dot(B1_Transpose, C_2)).dot(B3)
s32 = (np.dot(B2_Transpose, C_2)).dot(B3)
s33 = (np.dot(B3_Transpose, C_2)).dot(B3)
s34 = (np.dot(B4_Transpose, C_2)).dot(B3)
s41 = (np.dot(B1_Transpose, C_2)).dot(B4)
s42 = (np.dot(B2_Transpose, C_2)).dot(B4)
s43 = (np.dot(B3_Transpose, C_2)).dot(B4)
s44 = (np.dot(B4_Transpose, C_2)).dot(B4)

# print(f'\nK11 formula: {s11}')
sub_s11 = np.empty_like(s11)
sub_s12 = np.empty_like(s12)
sub_s13 = np.empty_like(s13)
sub_s14 = np.empty_like(s14)
sub_s21 = np.empty_like(s21)
sub_s22 = np.empty_like(s22)
sub_s23 = np.empty_like(s23)
sub_s24 = np.empty_like(s24)
sub_s31 = np.empty_like(s31)
sub_s32 = np.empty_like(s32)
sub_s33 = np.empty_like(s33)
sub_s34 = np.empty_like(s34)
sub_s41 = np.empty_like(s41)
sub_s42 = np.empty_like(s42)
sub_s43 = np.empty_like(s43)
sub_s44 = np.empty_like(s44)

'''
for i in range(s11.shape[0]):
    for j in range(s11.shape[1]):
        sub_s11[i, j] = s11[i, j].subs({inv_area: 1 / area, x_l: cor[0][0], x_h: cor[1][0], y_l: cor[0][1], y_h: cor[2][1], v: v_value})
'''
# K with symbols only
K = [
    s11, s21, s31, s41,
    s12, s22, s32, s42,
    s13, s23, s33, s43,
    s14, s24, s34, s44
]

# K with values substituted
K_subbed = [
    sub_s11, sub_s21, sub_s31, sub_s41,
    sub_s12, sub_s22, sub_s32, sub_s42,
    sub_s13, sub_s23, sub_s33, sub_s43,
    sub_s14, sub_s24, sub_s34, sub_s44
]

# Loop through K and K_subbed in parallel, substitute the values, and insert into K_subbed
for i1, i2 in zip(K, K_subbed):
# Loop through s11 and substitute the values
    for j in range(i1.shape[0]):
        for k in range(i1.shape[1]):
            # sub_s11[j, k] = s11[j , k].subs({inv_area: 1/area, x_l: cor[0][0], x_h: cor[1][0], y_l: cor[0][1], y_h: cor[2][1], v: v_value})
            i2[j, k] = i1[j, k].subs({inv_area: 1/area, x_l: cor[0][0], x_h: cor[1][0], y_l: cor[0][1], y_h: cor[2][1], v: v_value})


print(f'\nK11 substituted formula: \n{C_1} x \n{sub_s11}')
print(f'\nK21 substituted formula: \n{C_1} x \n{sub_s21}')

# Substitute E and Loop through K_subbed and multiply by C_1
C_1 = C_1.subs({E: E_value, v: v_value})
print(f'\nC_1: \n{C_1}')
K_subbed[0] = np.dot(C_1, sub_s11)
K_subbed[1] = np.dot(C_1, sub_s21)
K_subbed[2] = np.dot(C_1, sub_s31)
K_subbed[3] = np.dot(C_1, sub_s41)
K_subbed[4] = np.dot(C_1, sub_s12)
K_subbed[5] = np.dot(C_1, sub_s22)
K_subbed[6] = np.dot(C_1, sub_s32)
K_subbed[7] = np.dot(C_1, sub_s42)
K_subbed[8] = np.dot(C_1, sub_s13)
K_subbed[9] = np.dot(C_1, sub_s23)
K_subbed[10] = np.dot(C_1, sub_s33)
K_subbed[11] = np.dot(C_1, sub_s43)
K_subbed[12] = np.dot(C_1, sub_s14)
K_subbed[13] = np.dot(C_1, sub_s24)
K_subbed[14] = np.dot(C_1, sub_s34)
K_subbed[15] = np.dot(C_1, sub_s44)
int_s11 = np.empty_like(sub_s11)
int_s21 = np.empty_like(sub_s21)
int_s31 = np.empty_like(sub_s31)
int_s41 = np.empty_like(sub_s41)
int_s12 = np.empty_like(sub_s12)
int_s22 = np.empty_like(sub_s22)
int_s32 = np.empty_like(sub_s32)
int_s42 = np.empty_like(sub_s42)
int_s13 = np.empty_like(sub_s13)
int_s23 = np.empty_like(sub_s23)
int_s33 = np.empty_like(sub_s33)
int_s43 = np.empty_like(sub_s43)
int_s14 = np.empty_like(sub_s14)
int_s24 = np.empty_like(sub_s24)
int_s34 = np.empty_like(sub_s34)
int_s44 = np.empty_like(sub_s44)

K_integrated = [
    int_s11, int_s21, int_s31, int_s41,
    int_s12, int_s22, int_s32, int_s42,
    int_s13, int_s23, int_s33, int_s43,
    int_s14, int_s24, int_s34, int_s44
]

print(f'\nK11 Final: \n{sub_s11}')
print(f'\nK21 Final: \n{sub_s21}')

# intergrals of the stiffness matrix
'''
for j in range(sub_s21.shape[0]):
    for k in range(sub_s21.shape[1]):
        sub_s21[j, k] = integrate(sub_s21[j, k], (x, cor[0][0], cor[1][0]), (y, cor[0][1], cor[2][1]))

print(f'\nK21 integrated: \n{sub_s21.transpose()}')
'''

for i1, i2 in zip(K_subbed, K_integrated):
    for j in range(i1.shape[0]):
        for k in range(i1.shape[1]):
            i2[j, k] = integrate(i1[j, k], (x, cor[0][0], cor[1][0]), (y, cor[0][1], cor[2][1]))
# sub_s21 = [integrate(i, (x, cor[0][0], cor[1][0]), (y, cor[0][1], cor[2][1])) for i in sub_s21]
print(f'\nK11 integrated: \n{int_s11}')
print(f'\nK21 integrated: \n{int_s21}')
print(f'\nK31 integrated: \n{int_s31}')
print(f'\nK41 integrated: \n{int_s41}')
print(f'\nK12 integrated: \n{int_s12}')
print(f'\nK22 integrated: \n{int_s22}')
print(f'\nK32 integrated: \n{int_s32}')
print(f'\nK42 integrated: \n{int_s42}')
print(f'\nK13 integrated: \n{int_s13}')
print(f'\nK23 integrated: \n{int_s23}')
print(f'\nK33 integrated: \n{int_s33}')
print(f'\nK43 integrated: \n{int_s43}')
print(f'\nK14 integrated: \n{int_s14}')
print(f'\nK24 integrated: \n{int_s24}')
print(f'\nK34 integrated: \n{int_s34}')
print(f'\nK44 integrated: \n{int_s44}')





# test_list = [1,2,3,4]
# test_list = test_list * 2
# print(f'\nTest_list: {test_list}')
# test_list = [1,2,3,4]
# for i in test_list:
#     i = i * 2
# print(f'\nTest_list_looped: {test_list}')
# test_list = [1,2,3,4]
# test_list_output = []
# for i in test_list:
#     i = np.dot(i, 2)
#     test_list_output.append(i)
# print(f'\nTest_list_looped_np: {test_list}')
# print(f'\nTest_list_looped_np_output: {test_list_output}')
# test_array = np.array([1,2,3,4])
# test_array = test_array * 2
# print(f'\nTest_array: {test_array}')
