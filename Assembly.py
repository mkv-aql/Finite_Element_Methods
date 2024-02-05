__author__ = 'mkv-aql'
import numpy as np
from fractions import Fraction
import sympy as sp
import Function_Module as sm

nodes = 4
displacement = nodes * 4
boundaries_dataframe = ['3']
boundaries_numpy = [3]
print(boundaries_numpy[0])

#Create Stiffness matrix with 2 different angle nodes 1 and 2
K1G = sm.Stiffness_Matrix_1(45, 0)
# print(K1G)

#Create Stiffness matrix with single angle or both same angle
K3G = sm.Stiffness_Matrix_2(45)
# print(K2G + K3G)

#Remove boundary conditions numpy
print(sm.Boundary_Conditions_np(K3G, boundaries_numpy))


#Inserting variables of displacement, will turn into DataFrame type
K3G_new = sm.Insert_Variables(K3G)
print(K3G_new)
K3G_new_2 = sm.Insert_Variables(K3G, variable=['a','b','c','d'])
print(K3G_new_2)

#Remove boundary conditions DataFrame
print(sm.Boundary_Conditions_df(K3G_new, boundaries_dataframe))



#Assemble matrix
boundaries_dataframe_1 = ['u1', 'v1', 'u3', 'v3', 'v4']

A = sm.Stiffness_Matrix_2(-30, k = 20)
B = sm.Stiffness_Matrix_2(45, k = 15)
C = sm.Stiffness_Matrix_2(0, k = 10)
# print(A,B,C)

A = sm.Insert_Variables(A, variable = ['u1','v1','u2','v2'])
B = sm.Insert_Variables(B, variable = ['u2','v2','u3','v3'])
C = sm.Insert_Variables(C, variable = ['u2','v2','u4','v4'])
# print(A,B,C)

A = sm.Boundary_Conditions_df(A, boundaries_dataframe_1)
B = sm.Boundary_Conditions_df(B, boundaries_dataframe_1)
C = sm.Boundary_Conditions_df(C, boundaries_dataframe_1)

print(sm.Combine_Matrix_df([A,B,C]))


#Assemble matrix 2
k = 20
boundaries_dataframe_2 = ['u1', 'v1', 'v3']
K1= sm.Stiffness_Matrix_2(90, k = k)
K2 = sm.Stiffness_Matrix_2(-90, k = k)
K3 = sm.Stiffness_Matrix_2(150, k = k)
K4 = sm.Stiffness_Matrix_2(-150, k = k)
K5 = sm.Stiffness_Matrix_2(150, k = k)
K6 = sm.Stiffness_Matrix_2(-150, k = k)
K7 = sm.Stiffness_Matrix_2(90, k = k)

K1 = sm.Insert_Variables(K1, variable = ['u1','v1','u2','v2'])
K2 = sm.Insert_Variables(K2, variable = ['u2','v2','u3','v3'])
K3 = sm.Insert_Variables(K3, variable = ['u1','v1','u4','v4'])
K4 = sm.Insert_Variables(K4, variable = ['u2','v2','u4','v4'])
K5 = sm.Insert_Variables(K5, variable = ['u2','v2','u5','v5'])
K6 = sm.Insert_Variables(K6, variable = ['u3','v3','u5','v5'])
K7 = sm.Insert_Variables(K7, variable = ['u4','v4','u5','v5'])

K1 = sm.Boundary_Conditions_df(K1, boundaries_dataframe_2)
K2 = sm.Boundary_Conditions_df(K2, boundaries_dataframe_2)
K3 = sm.Boundary_Conditions_df(K3, boundaries_dataframe_2)
K4 = sm.Boundary_Conditions_df(K4, boundaries_dataframe_2)
K5 = sm.Boundary_Conditions_df(K5, boundaries_dataframe_2)
K6 = sm.Boundary_Conditions_df(K6, boundaries_dataframe_2)
K7 = sm.Boundary_Conditions_df(K7, boundaries_dataframe_2)

print(sm.Combine_Matrix_df([K1,K2,K3,K4,K5,K6,K7]))