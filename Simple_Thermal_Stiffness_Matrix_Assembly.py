__author__ = 'mkv-aql'
import math
import numpy as np
import pandas as pd
from sympy import Symbol, cos, factor, symbols
import copy # to create deep copy of the list, like a quick instance without using classes and objects

# List of local nodes, bottom left, bottom right, top right, top left
local_nodes = [1, 2, 3, 4] # Parameter
local_nodes_matrix = np.zeros((4,4)) # Parameter
Assembed_K = np.zeros((8,8)) # Parameter
# Stiffness_K_Value = 1 # Parameter for example 1
Stiffness_K_Value = 500.0556731 # Parameter for example 2

# List of elements, 1st element, 2nd element, 3rd element, etc. / create 3 elements instance of local nodes
elements = [copy.deepcopy(local_nodes) for i in range(3)] # Parameter
e = elements

# List of global nodes, bottom left, bottom right, top right, top left
global_nodes = [1, 2, 3, 4, 5, 6, 7, 8] # Parameter
gn = global_nodes

print(elements[1][3])

# Example 1 of stiffness matrix
K = np.array([[0.7556, 0.0778, -0.3778, -0.4556],
              [0.0778, 0.7556, -0.4556, -0.3778],
              [-0.3778, -0.4556, 0.7556, 0.0778],
              [-0.4556, -0.3778, 0.0778, 0.7556]])

# Example 2 of stiffness matrix
K = np.array([[0.8981, 0.2269, -0.4491, -0.6759],
             [0.2269, 0.8981, -0.6759, -0.4491],
             [-0.4491, -0.6759, 0.8981, 0.2269],
             [-0.6759, -0.4491, 0.2269, 0.8981]])

# Assign the stiffness matrix to 1 element
local_nodes_matrix[:] = K
print(f'Stiffness matrix for 1 element: \n{local_nodes_matrix}')

# Assign the stiffness matrix to all elements
local_nodes_matrix_1 = np.zeros((4,4))
local_nodes_matrix_2 = np.zeros((4,4))
local_nodes_matrix_3 = np.zeros((4,4))
local_nodes_matrix_1[:] = K
local_nodes_matrix_2[:] = K
local_nodes_matrix_3[:] = K
K_1 = local_nodes_matrix_1
K_2 = local_nodes_matrix_2
K_3 = local_nodes_matrix_3


# Example of Thermal U matrix
u1, u2, u3, u4, u5, u6, u7, u8 = symbols('u1, u2, u3, u4, u5, u6, u7, u8')
U = np.array([u1, u2, u3, u4, u5, u6, u7, u8])

# Local Thermal U
Local_U_1 = ['u1', 'u2', 'u3', 'u4']
Local_U_2 = ['u3', 'u4', 'u6', 'u7']
Local_U_3 = ['u4', 'u5', 'u7', 'u8']
Local_U_1_g = ['u1', 'u2', 'u4', 'u3']
Local_U_2_g = ['u3', 'u4', 'u7', 'u6']
Local_U_3_g = ['u4', 'u5', 'u8', 'u7']

Local_Us = [Local_U_1, Local_U_2, Local_U_3]
Local_Us_g = [Local_U_1_g, Local_U_2_g, Local_U_3_g]

'''
# Basic of rearranging the global U 1 by 1
# Create a mapping from the initial order to the new order
index_map_1 = [Local_U_1.index(u) for u in Local_U_1_g]
# Add 1 to make it human readable
index_map_1_show = [i+1 for i in index_map_1] # ONly for human readable, not used in calculation
print(f'\nIndex map 1: \n{index_map_1_show}')

# test = ['u1', 'u2', 'u3', 'u4', 'u5', 'u6', 'u7', 'u8']
# index_map = [Local_U_1.index(u) for u in test]

# Rearrange the matrix K according to the new order of U global
K_1 = K[:, index_map_1][index_map_1, :]
print(f'\nNew K_1: \n{K_1}')
'''
index_map = [[],[],[]]
symbol_map = [[],[],[]]
n = 0
for Local_U, Local_U_g in zip(Local_Us, Local_Us_g): #Zip the Lists: Use zip(Local_Us, Local_Us_g) to iterate over corresponding pairs of Local_U and Local_U_g together. This ensures you're mapping elements of the correct corresponding lists.
    for u in Local_U_g:
        if u in Local_U:
            index_map[n].append(Local_U.index(u)) # Will get the index/indices values, +1 to make it human readable
            symbol_map[n].append(u) # Append into the symbol_map per dimension
        else:
            index_map[n].append(None)
            symbol_map[n].append(None)
    n += 1

print(f'\nMapping Connectivity of Element Nodes: \n{np.array(index_map)}') #np.array to get a nicer matrix view
print(f'\nMapping Connectivity of Elements: \n{np.array(symbol_map)}')


print(f'\nThermal U: \n{U}')

# Rearrange the matrix K according to the new order of U global
K_1 = K[:, index_map[0]][index_map[0], :]
print(f'\nNew K_1: \n{K_1}')
K_2 = K[:, index_map[1]][index_map[1], :]
print(f'\nNew K_2: \n{K_2}')
K_3 = K[:, index_map[2]][index_map[2], :]
print(f'\nNew K_3: \n{K_3}')

def add_to_global_matrix(Assembed_K, local_matrix, connections):
    for i in range(4):  # Loop through rows of the 4x4 matrix
        for j in range(4):  # Loop through columns of the 4x4 matrix
            Assembed_K[connections[i]][connections[j]] += local_matrix[i][j]

# Create a mapping from symbol to number
symbol_to_number = {'u1': 0, 'u2': 1, 'u3': 2, 'u4': 3, 'u5': 4, 'u6': 5, 'u7': 6, 'u8': 7}

# Convert the symbols array to a numerical array
number_array = np.vectorize(symbol_to_number.get)(symbol_map)
print(f'\nNumber array element based: \n{number_array}')

number_array = np.array([[0,1,2,3],[2,3,5,6],[3,4,6,7]])
print(f'\nNumber array global based: \n{number_array}')
add_to_global_matrix(Assembed_K, K_1, number_array[0])
add_to_global_matrix(Assembed_K, K_2, number_array[1])
add_to_global_matrix(Assembed_K, K_3, number_array[2])

print(f'\nAssembed K: \n{Assembed_K}')

# Multiply the stiffness value to the global matrix
Assembed_K *= Stiffness_K_Value

for i in range(8):
    for j in range(8):
        Assembed_K[i][j] = round(Assembed_K[i][j], 2)

print(f'\nAssembed K multiplied by stiffness value: \n{Assembed_K}')