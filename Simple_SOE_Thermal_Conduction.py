__author__ = 'mkv-aql'
import numpy as np
from sympy import symbols, Eq, solve
'''
[[ 449.1   113.46 -337.99 -224.58    0.      0.      0.      0.  ]
 [ 113.46  449.1  -224.58 -337.99    0.      0.      0.      0.  ]
 [-337.99 -224.58  898.2   226.93    0.   -337.99 -224.58    0.  ]
 [-224.58 -337.99  226.93 1347.3   113.46 -224.58 -675.98 -224.58]
 [   0.      0.      0.    113.46  449.1     0.   -224.58 -337.99]
 [   0.      0.   -337.99 -224.58    0.    449.1   113.46    0.  ]
 [   0.      0.   -224.58 -675.98 -224.58  113.46  898.2   113.46]
 [   0.      0.      0.   -224.58 -337.99    0.    113.46  449.1 ]]
'''
# Example 2
K = np.array([[449.1, 113.46, -337.99, -224.58, 0, 0, 0, 0],
              [113.46, 449.1, -224.58, -337.99, 0, 0, 0, 0],
              [-337.99, -224.58, 898.2, 226.93, 0, -337.99, -224.58, 0],
              [-224.58, -337.99, 226.93, 1347.3, 113.46, -224.58, -675.98, -224.58],
              [0, 0, 0, 113.46, 449.1, 0, -224.58, -337.99],
              [0, 0, -337.99, -224.58, 0, 449.1, 113.46, 0],
              [0, 0, -224.58, -675.98, -224.58, 113.46, 898.2, 113.46],
              [0, 0, 0, -224.58, -337.99, 0, 113.46, 449.1]])

'''
[[ 0.76  0.08 -0.46 -0.38  0.    0.    0.    0.  ]
 [ 0.08  0.76 -0.38 -0.46  0.    0.    0.    0.  ]
 [-0.46 -0.38  1.51  0.16  0.   -0.46 -0.38  0.  ]
 [-0.38 -0.46  0.16  2.27  0.08 -0.38 -0.91 -0.38]
 [ 0.    0.    0.    0.08  0.76  0.   -0.38 -0.46]
 [ 0.    0.   -0.46 -0.38  0.    0.76  0.08  0.  ]
 [ 0.    0.   -0.38 -0.91 -0.38  0.08  1.51  0.08]
 [ 0.    0.    0.   -0.38 -0.46  0.    0.08  0.76]]
 '''
#Example 1, uncomment to use
# K = np.array([[0.76, 0.08, -0.46, -0.38, 0, 0, 0, 0],
#                 [0.08, 0.76, -0.38, -0.46, 0, 0, 0, 0],
#                 [-0.46, -0.38, 1.51, 0.16, 0, -0.46, -0.38, 0],
#                 [-0.38, -0.46, 0.16, 2.27, 0.08, -0.38, -0.91, -0.38],
#                 [0, 0, 0, 0.08, 0.76, 0, -0.38, -0.46],
#                 [0, 0, -0.46, -0.38, 0, 0.76, 0.08, 0],
#                 [0, 0, -0.38, -0.91, -0.38, 0.08, 1.51, 0.08],
#                 [0, 0, 0, -0.38, -0.46, 0, 0.08, 0.76]])

K_hat = K #copies
K_hat_F = K #copies

F = np.zeros(8)

u1, u2, u3, u4, u5, u6, u7, u8 = symbols('u1, u2, u3, u4, u5, u6, u7, u8') #parameter
U = np.array([u1, u2, u3, u4, u5, u6, u7, u8]) # parameter
known_values = {u8: 20} #parameter
boundary_conditions = {u1: 0, u2: 0} #parameter
unknown_U = np.array([u for u in U if u not in known_values.keys() and u not in boundary_conditions.keys()]) # [u3 u4 u5 u6 u7]
# print(known_values.keys()) #dict_keys([u8])
# # Get index value of u8 using index()
# u8 = U.tolist().index(u8)
# print(f'\nIndex of u8: {u8}') # 7
# # Get index value of u8 from known_values in U
# u8 = U.tolist().index(list(known_values.keys())[0])
# print(f'\nIndex of u8: {u8}') #7
# Get value of u8 from known_values
# u8 = list(known_values.values())[0]
# print(f'\nValue of u8: {u8}') #20
# u1 = U.tolist().index(list(boundary_conditions.keys())[0])
# print(f'\nIndex of u1: {u1}')
# u2 = U.tolist().index(list(boundary_conditions.keys())[1])
# print(f'\nIndex of u1: {u1}')
# len of boundary_conditions
# print(f'\nLength of boundary_conditions: {len(boundary_conditions)}')

print(f'\n K matrix Original: \n{K}')

#Reduce K and F from boundary conditions and known values
for i in range(len(known_values)):
    K_hat = np.delete(K_hat, U.tolist().index(list(known_values.keys())[0]), axis=0) #Cannot use i, must use 0, because after the first deletion, when i = 1, it will delete the 2nd row, not the 1st row
    K_hat = np.delete(K_hat, U.tolist().index(list(known_values.keys())[0]), axis=1) #Cannot use i, must use 0, because after the first deletion, when i = 1, it will delete the 2nd row, not the 1st row
    F = np.delete(F, U.tolist().index(list(known_values.keys())[i]), axis=0)
# K = np.delete(K, U.tolist().index(list(known_values.keys())[0]), axis=0)
# K = np.delete(K, U.tolist().index(list(known_values.keys())[0]), axis=1)
for i in range(len(boundary_conditions)):
    K_hat = np.delete(K_hat, U.tolist().index(list(boundary_conditions.keys())[0]), axis=0) #Cannot use i, must use 0, because after the first deletion, when i = 1, it will delete the 2nd row, not the 1st row
    K_hat = np.delete(K_hat, U.tolist().index(list(boundary_conditions.keys())[0]), axis=1) #Cannot use i, must use 0, because after the first deletion, when i = 1, it will delete the 2nd row, not the 1st row
    F = np.delete(F, U.tolist().index(list(boundary_conditions.keys())[i]), axis=0)
# K = np.delete(K, U.tolist().index(list(boundary_conditions.keys())[0]), axis=0)
# K = np.delete(K, U.tolist().index(list(boundary_conditions.keys())[0]), axis=1)
# K = np.delete(K, U.tolist().index(list(boundary_conditions.keys())[1]), axis=0)
# K = np.delete(K, U.tolist().index(list(boundary_conditions.keys())[1]), axis=1)
print(f'\nReduced K matrix: \n{K_hat}')
print(f'\nReduced F: \n{F}')

# Get the known values of U, multiply to its corresponding row in K and subtract from F
    #Create the new K_hat_F, the known values of U in K, with reduced boundary conditions and known values
for i in range(len(known_values)):
    K_hat_F = np.delete(K_hat_F[U.tolist().index(list(known_values.keys())[i])], U.tolist().index(list(known_values.keys())[i]), axis=0) #[   0.      0.      0.   -224.58 -337.99    0.    113.46]
for i in range(len(boundary_conditions)):
    K_hat_F = np.delete(K_hat_F, U.tolist().index(list(boundary_conditions.keys())[i]), axis=0) #[   0.   -224.58 -337.99    0.    113.46]
print(f'\nReduced K for F: \n{K_hat_F}')
# print(f'\nTest: \n{np.dot(K_hat_F, (list(known_values.values())[0]))}') #[    0.   4491.6  6759.8     0.  -2269.2]
F = F - np.dot(K_hat_F, list(known_values.values())[0])
print(f'\nNew F: \n{F}')


# Solve for the unknown values of U
U = np.linalg.solve(K_hat, F)
print(f'\nU: \n{U}')
np.allclose(np.dot(K_hat, U), F) #Check if the solution is correct

#Dataframes
import pandas as pd

# df = pd.DataFrame(K, columns=['u1', 'u2', 'u3', 'u4', 'u5', 'u6', 'u7', 'u8'])
df = pd.DataFrame(K_hat, columns=[unknown_U])
df2 = pd.DataFrame(F)
# show column names
print(df.columns)

df = df.to_numpy()
df2 = df2.to_numpy()

print(f'\nU: \n{np.linalg.solve(df, df2)}')

