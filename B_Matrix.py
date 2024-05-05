__author__ = 'mkv-aql'

import numpy as np
import Function_Module as sm
import Basic_Values as bv

Node_Coordinates = [[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]]

'''
com, dh_drs, jac, _, inv_jac = sm.Jacobian_Gauss_Legendre(Node_Coordinates)

print(com[2]) #Gauss point 3
print(dh_drs[2]) #Derivative of shape function for 3rd Gauss point
print(jac[2]) #Jacobian for 3rd Gauss point
print(inv_jac[2]) #Inverse Jacobian for 3rd Gauss point

# Parameters
chosen_node = 2 # B matrix of node 2
chosen_gauss_point = 3 # B matrix of 3rd Gauss point

#Formula
result = np.dot(inv_jac[chosen_gauss_point - 1], dh_drs[chosen_gauss_point - 1][chosen_node - 1]) # formula
dh_dx, dh_dy = result[0], result[1]

# Insert into B matrix
b_matrix = [[dh_dx, 0],[0, dh_dy],[dh_dy, dh_dx]]

print("B Matrix = ", b_matrix)
'''

def B_Matrix(node_coordinates, integration_point, chosen_node):
    com, dh_drs, jac, _, inv_jac = sm.Jacobian_Gauss_Legendre(node_coordinates, integration_point)
    result = np.dot(inv_jac[integration_point - 1], dh_drs[integration_point - 1][chosen_node - 1]) # formula
    dh_dx, dh_dy = result[0], result[1]
    b_matrix = [[dh_dx, 0],[0, dh_dy],[dh_dy, dh_dx]]
    return b_matrix

print(B_Matrix(Node_Coordinates, 3, 1)) # B matrix of node 1
print(B_Matrix(Node_Coordinates, 3, 2)) # B matrix of node 2
print(B_Matrix(Node_Coordinates, 3, 3)) # B matrix of node 3
print(B_Matrix(Node_Coordinates, 3, 4)) # B matrix of node 4

# convert to dataframes
B_Matrix_df = sm.Insert_Variables(np.array(B_Matrix(Node_Coordinates, 3, 1)), columns = ['B1', 'B1'], index = ['1', '2', '3'])
print(B_Matrix_df)
B_Matrix_df = sm.Insert_Variables(np.array(B_Matrix(Node_Coordinates, 3, 2)), columns = ['B2', 'B2'], index = ['1', '2', '3'])
print(B_Matrix_df)
B_Matrix_df = sm.Insert_Variables(np.array(B_Matrix(Node_Coordinates, 3, 3)), columns = ['B3', 'B3'], index = ['1', '2', '3'])
print(B_Matrix_df)
B_Matrix_df = sm.Insert_Variables(np.array(B_Matrix(Node_Coordinates, 3, 4)), columns = ['B4', 'B4'], index = ['1', '2', '3'])
print(B_Matrix_df)

# B matrix with only chosen integration points n (to get rg and sg location of gauss point) and node coordinates
def B_Matrix_Assembly(node_coordinates, integration_point):
    com = bv.gauss_quadrature_r_s_combination[integration_point]
    rg_sg_location = 1

    _, _, jac, det_jac, inv_jac = sm.Jacobian_Gauss_Legendre(node_coordinates, integration_point)

    return 0

def B_Matrix_Assembly_Of_Combination(node_coordinates, integration_point, combination):
    com, dh_drs, jac, _, inv_jac = sm.Jacobian_Gauss_Legendre(node_coordinates, integration_point)
    final_b_matrix = []
    for node in node_coordinates:
        #index = node_coordinates.index(node)
        b_matrix = B_Matrix(node_coordinates, integration_point, node)
        final_b_matrix.append(b_matrix)


    return final_b_matrix

print(B_Matrix_Assembly_Of_Combination(Node_Coordinates, integration_point = 2, combination = 1))


