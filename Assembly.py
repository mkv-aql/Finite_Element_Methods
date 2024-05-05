__author__ = 'mkv-aql'
import numpy as np
from fractions import Fraction
import sympy as sp
import Function_Module as sm

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
Initial_DOF_1 = ['u1', 'v1', 'u2', 'v2', 'u3', 'v3', 'u4', 'v4'] #Initial DOF of system
Force_Vector_1 = [0, 0, 10] #Force Vector, adjust size to match DOF

DOF_1 = sm.Final_DOF(Initial_DOF_1, boundaries_dataframe_1)

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

Assembled_Matrix = sm.Combine_Matrix_df([A,B,C])
Assembled_Matrix = sm.Reorder_Matrix(Assembled_Matrix, DOF_1)

print("\nAssembled Matrix_1")
print(Assembled_Matrix)
print(sm.System_Of_Equation(Assembled_Matrix, Force_Vector_1))

#Assemble matrix 2
k = 20
boundaries_dataframe_2 = ['u1', 'v1', 'v3']
Initial_DOF_2 = ['u1', 'v1', 'u2', 'v2', 'u3', 'v3', 'u4', 'v4', 'u5', 'v5'] #Initial DOF of system
Force_Vector_2 = [0, -1, 0, 0, 0, 0, 0]

DOF_2 = sm.Final_DOF(Initial_DOF_2, boundaries_dataframe_2)
# print(DOF_2)

K1 = sm.Stiffness_Matrix_2(0, k = k)
K2 = sm.Stiffness_Matrix_2(0, k = k)
K3 = sm.Stiffness_Matrix_2(60, k = k)
K4 = sm.Stiffness_Matrix_2(120, k = k)
K5 = sm.Stiffness_Matrix_2(60, k = k)
K6 = sm.Stiffness_Matrix_2(120, k = k)
K7 = sm.Stiffness_Matrix_2(0, k = k)

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

Assembled_Matrix_2 = sm.Combine_Matrix_df([K1,K2,K3,K4,K5,K6,K7])
Assembled_Matrix_2 = sm.Reorder_Matrix(Assembled_Matrix_2, DOF_2)
# print(Assembled_Matrix_2.to_numpy())

print("\nAssembled Matrix_2")
print(Assembled_Matrix_2)
print(sm.System_Of_Equation(Assembled_Matrix_2, Force_Vector_2))

#Assemble matrix 3
k_3 = 20
f_1 = -1
f_2 = 0.5
boundaries_dataframe_3 = ['u5', 'v5', 'u6', 'v6']
Initial_DOF_3 = ['u1', 'v1', 'u2', 'v2', 'u3', 'v3', 'u4', 'v4', 'u5', 'v5', 'u6', 'v6'] #Initial DOF of system
#
Force_Vector_3 = [0, f_1, f_2, 0, 0, 0, 0, 0] # Force Vector, adjusted size to match DOF

DOF_3 = sm.Final_DOF(Initial_DOF_3, boundaries_dataframe_3)

#Ascending order of node, lowest node will be pivot, but if already horizontal, then just cos 180 or sin 180 with displacement label will be ascending order
# (lowest node must be on the left side of diagram)
I = sm.Stiffness_Matrix_2(180, k = k_3)
II = sm.Stiffness_Matrix_2(180, k = k_3)
III = sm.Stiffness_Matrix_2(180, k = k_3)
IV = sm.Stiffness_Matrix_2(180, k = k_3)
V = sm.Stiffness_Matrix_2(-90, k = k_3) #Vertical element, lowest node will be pivot
VI = sm.Stiffness_Matrix_2(-90, k = k_3) #Vertical element, lowest node will be pivot
VII = sm.Stiffness_Matrix_2(135, k = k_3)
VIII = sm.Stiffness_Matrix_2(-135, k = k_3)
IX = sm.Stiffness_Matrix_2(135, k = k_3)
X = sm.Stiffness_Matrix_2(-135, k = k_3)

I = sm.Insert_Variables(I, variable = ['u3','v3','u5','v5'])
II = sm.Insert_Variables(II, variable = ['u1','v1','u3','v3'])
III = sm.Insert_Variables(III, variable = ['u4','v4','u6','v6'])
IV = sm.Insert_Variables(IV, variable = ['u2','v2','u4','v4'])
V = sm.Insert_Variables(V, variable = ['u3','v3','u4','v4'])
VI = sm.Insert_Variables(VI, variable = ['u1','v1','u2','v2'])
VII = sm.Insert_Variables(VII, variable = ['u4','v4','u5','v5'])
VIII = sm.Insert_Variables(VIII, variable = ['u3','v3','u6','v6'])
IX = sm.Insert_Variables(IX, variable = ['u2','v2','u3','v3'])
X = sm.Insert_Variables(X, variable = ['u1','v1','u4','v4'])

I = sm.Boundary_Conditions_df(I, boundaries_dataframe_3)
II = sm.Boundary_Conditions_df(II, boundaries_dataframe_3)
III = sm.Boundary_Conditions_df(III, boundaries_dataframe_3)
IV = sm.Boundary_Conditions_df(IV, boundaries_dataframe_3)
V = sm.Boundary_Conditions_df(V, boundaries_dataframe_3)
VI = sm.Boundary_Conditions_df(VI, boundaries_dataframe_3)
VII = sm.Boundary_Conditions_df(VII, boundaries_dataframe_3)
VIII = sm.Boundary_Conditions_df(VIII, boundaries_dataframe_3)
IX = sm.Boundary_Conditions_df(IX, boundaries_dataframe_3)
X = sm.Boundary_Conditions_df(X, boundaries_dataframe_3)

Assembled_Matrix_3 = sm.Combine_Matrix_df([I,II,III,IV,V,VI,VII,VIII,IX,X])
Assembled_Matrix_3 = sm.Reorder_Matrix(Assembled_Matrix_3, DOF_3)
print(" \nAssembled Matrix_3")
print(Assembled_Matrix_3)
print(sm.System_Of_Equation(Assembled_Matrix_3, Force_Vector_3))