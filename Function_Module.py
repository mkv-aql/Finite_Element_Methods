__author__ = 'mkv-aql'
import numpy as np
from fractions import Fraction
import sympy as sp
import pandas as pd
from sympy import Symbol, cos, factor, symbols, simplify
import math
import Basic_Values as bv

#Print whole matrix if exceed 5x5
pd.set_option('display.max_columns', None)

def Stiffness_Matrix_1(alpha, beta, k = 1):
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
    return stiffness_matrix * k

def Stiffness_Matrix_2(angle, k = 1):
    angle = angle / 360 * 2 * np.pi

    K_Matrix = np.array([[1, -1],
                         [-1, 1]])

    Transformation = np.array([[np.cos(angle), np.sin(angle), 0, 0],
                               [0, 0, np.cos(angle), np.sin(angle)]])

    Transformation_transpose = Transformation.transpose()

    Stiffness_Matrix_Transformed = Transformation_transpose.dot(K_Matrix.dot(Transformation))

    return Stiffness_Matrix_Transformed * k

def Transformation(angle):
    Transformation = np.array([[np.cos(angle), np.sin(angle), 0, 0],
                               [0, 0, np.cos(angle), np.sin(angle)]])

    Transformation_transpose = Transformation.transpose()

    return Transformation, Transformation_transpose

def Boundary_Conditions_np(matrix, boundaries):
    for bound in boundaries:
        matrix_1 = np.delete(matrix, bound, axis=0) #Remove row
        matrix_2 = np.delete(matrix_1, bound, axis=1) #Remove column

    return matrix_2

def Insert_Variables(matrix, variable = ['1', '2', '3', '4']):
    df = pd.DataFrame(matrix, columns=variable, index=variable)
    return df

def Insert_Variables(matrix, columns, index):
    df = pd.DataFrame(matrix, columns=columns, index=index)
    return df

def Boundary_Conditions_df(dataframe, boundaries):
    for bound in boundaries:
        try:
            dataframe = dataframe.drop(columns=bound, index = bound)
        except:
            pass

    return dataframe

def Combine_Matrix_df(matrix_list):
    # Initialize an empty DataFrame to combine everything
    combined_df = pd.DataFrame()
    # List of all DataFrames for iteration
    dfs = matrix_list

    for df in dfs:
        # If combined_df is empty, simply assign the first DataFrame
        if combined_df.empty:
            combined_df = df
        else:
            # For subsequent DataFrames, merge with combined_df on index (left_index and right_index)
            # using 'outer' join to ensure all variables are included
            combined_df = combined_df.merge(df, left_index=True, right_index=True, how='outer', suffixes=('', '_dup'))

            # After merging, columns with '_dup' are duplicates; add these to the original columns and drop '_dup' columns
            for col in df.columns:
                if col + '_dup' in combined_df.columns:
                    combined_df[col] = combined_df[col].fillna(0) + combined_df[col + '_dup'].fillna(0)
                    combined_df.drop(col + '_dup', axis=1, inplace=True)

    # Ensure the final DataFrame has no missing values (fill with 0 or appropriate value)
    combined_df.fillna(0, inplace=True)

    return combined_df

def Final_DOF(all_dof, boundaries):
    unique = set(all_dof) ^ set(boundaries)  # Remove boundary conditions from DOF
    Final_DOF = [x for x in all_dof if x in unique]  # Final DOF after removing boundary conditions

    return Final_DOF

def Reorder_Matrix(matrix, order):
    matrix = matrix.reindex(index=order, columns=order)
    return matrix

def System_Of_Equation(assembled_matrix, force_vector):
    if type(assembled_matrix) == 'pandas.core.frame.DataFrame':
        assembled_matrix = assembled_matrix.to_numpy()
    else:
        pass
    soe = np.linalg.solve(assembled_matrix, force_vector)

    return soe


def Jacobian_Numerical_Integration(coordinate):
    ##Jacobian calculation with Numerical Integration
    # Linear shape function for 4 node element in 2D (N1 = h1)
    r, s = symbols('r,s')
    h1 = (1 / 4) * (1 - r) * (1 - s)
    h2 = (1 / 4) * (1 + r) * (1 - s)
    h3 = (1 / 4) * (1 + r) * (1 + s)
    h4 = (1 / 4) * (1 - r) * (1 + s)

    # Node coordinates = [Node 1, Node 2, Node 3, Node 4]
    NC = coordinate

    x_rs = h1 * NC[0][0] + h2 * NC[1][0] + h3 * NC[2][0] + h4 * NC[3][0]
    y_rs = h1 * NC[0][1] + h2 * NC[1][1] + h3 * NC[2][1] + h4 * NC[3][1]

    print(x_rs)
    print(y_rs)

    print("dx = ", simplify(x_rs))
    print("dy = ", simplify(y_rs))

    # Jacobian, J = [[dx/dr, dy/dr], [dx/ds, dy/ds]]
    J = [[x_rs.diff(r), y_rs.diff(r)], [x_rs.diff(s), y_rs.diff(s)]]
    print("Jacobian = ", J)

    # Determinant of Jacobian
    J_det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    print("Det J = ", J_det)

    return J_det


def Jacobian_Gauss_Legendre(coordinate, integration_point):
    # Linear shape function for 4 node element in 2D (N1 = h1)
    r, s = symbols('r,s')
    h1 = (1 / 4) * (1 - r) * (1 - s)
    h2 = (1 / 4) * (1 + r) * (1 - s)
    h3 = (1 / 4) * (1 + r) * (1 + s)
    h4 = (1 / 4) * (1 - r) * (1 + s)

    # Node coordinates = [Node 1, Node 2, Node 3, Node 4]
    NC = coordinate

    # Parameters for Gauss Quadrature
    l_rg = math.sqrt(1 / 3)  # location of rg and sg
    weight = 1  # Wg
    c_r_s = [[-(l_rg), -(l_rg)], [(l_rg), -(l_rg)], [(l_rg), (l_rg)], [-(l_rg), (l_rg)]]  # Combination of r and s
    c_r_s = bv.gauss_quadrature_r_s_combination[integration_point] # Combination of r and s for 2nd Gauss point

    combination_of_r_s = []
    list_of_dh_drs = []
    output_jacobian = []
    output_jacobian_det = []
    output_jacobian_inv = []

    for r_s in c_r_s:
        #print("")
        # Iterate through the combination of r and s
        #print("combination of r and s = ", r_s)
        index = c_r_s.index(r_s)  # index of r_s in c_r_s

        # Bilinear shape function for 4 node element in 2D (N1 = h1), dh1/dr = h1_r, dh1/ds = h1_s, and so on.
        dh1_dr = h1.diff(r)
        dh1_ds = h1.diff(s)
        dh2_dr = h2.diff(r)
        dh2_ds = h2.diff(s)
        dh3_dr = h3.diff(r)
        dh3_ds = h3.diff(s)
        dh4_dr = h4.diff(r)
        dh4_ds = h4.diff(s)

        # substitute r and s values
        value_1, value_2 = c_r_s[index][0], c_r_s[index][1]
        dh1_dr = dh1_dr.subs({r: value_1, s: value_2})
        dh1_ds = dh1_ds.subs({r: value_1, s: value_2})
        dh2_dr = dh2_dr.subs({r: value_1, s: value_2})
        dh2_ds = dh2_ds.subs({r: value_1, s: value_2})
        dh3_dr = dh3_dr.subs({r: value_1, s: value_2})
        dh3_ds = dh3_ds.subs({r: value_1, s: value_2})
        dh4_dr = dh4_dr.subs({r: value_1, s: value_2})
        dh4_ds = dh4_ds.subs({r: value_1, s: value_2})
        dh_drs_list_current = [[dh1_dr, dh1_ds], [dh2_dr, dh2_ds], [dh3_dr, dh3_ds], [dh4_dr, dh4_ds]]

        # Jacobian, J = [[dx/dr, dy/dr], [dx/ds, dy/ds]] = [[A, B], [C, D]]
        dx_dr = dh1_dr * NC[0][0] + dh2_dr * NC[1][0] + dh3_dr * NC[2][0] + dh4_dr * NC[3][0]
        dx_ds = dh1_ds * NC[0][0] + dh2_ds * NC[1][0] + dh3_ds * NC[2][0] + dh4_ds * NC[3][0]
        dy_dr = dh1_dr * NC[0][1] + dh2_dr * NC[1][1] + dh3_dr * NC[2][1] + dh4_dr * NC[3][1]
        dy_ds = dh1_ds * NC[0][1] + dh2_ds * NC[1][1] + dh3_ds * NC[2][1] + dh4_ds * NC[3][1]

        Jac = [[dx_dr, dy_dr], [dx_ds, dy_ds]]
        #print("Jacobian = ", Jac)

        # Determinant of Jacobian
        Jac_det = dx_dr * dy_ds - dx_ds * dy_dr
        #print("Det J = ", Jac_det)

        # Convert to matrix
        Jac = np.array(Jac, dtype=float)

        # Inverse of Jacobian
        Jac_inv = np.linalg.inv(Jac)
        #print("Jacobian Inverse = ", Jac_inv)

        combination_of_r_s.append(r_s)
        list_of_dh_drs.append(dh_drs_list_current)
        output_jacobian.append(Jac)
        output_jacobian_det.append(Jac_det)
        output_jacobian_inv.append(Jac_inv)


    return combination_of_r_s, list_of_dh_drs, output_jacobian, output_jacobian_det, output_jacobian_inv

def Material_Matrix(E, v):
    c_matrix = np.array([[1, v, 0],
                     [v, 1, 0],
                     [0, 0, (1-v)/2]]) * E / (1 - v**2)
    c_matrix = np.round(c_matrix, 2)
    return c_matrix

def B_Matrix(Node_Coordinates, intgration_point, chosen_node):
    com, dh_drs, jac, _, inv_jac = Jacobian_Gauss_Legendre(Node_Coordinates, intgration_point)
    result = np.dot(inv_jac[intgration_point - 1], dh_drs[intgration_point - 1][chosen_node - 1]) # formula
    dh_dx, dh_dy = result[0], result[1]
    b_matrix = [[dh_dx, 0],[0, dh_dy],[dh_dy, dh_dx]]
    return b_matrix

a,_,_,_,_ = Jacobian_Gauss_Legendre([[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]], 2)
print("Com: ", a)
_,b,_,_,_ = Jacobian_Gauss_Legendre([[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]], 2)
print("dh_drs: ", b)
_,_,c,_,_ = Jacobian_Gauss_Legendre([[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]], 2)
print("Jac: ", c)
_,_,_,d,_ = Jacobian_Gauss_Legendre([[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]], 2)
print("Det Jac: ", d)
_,_,_,_,e = Jacobian_Gauss_Legendre([[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]], 2)
print("Jac Inv: ", e)

print(B_Matrix([[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]], 2, 2)) # B matrix of node 1

'''
def B_Matrix_Assembly(Node_Coordinates, integration_point):
    com, dh_drs, jac, _, inv_jac = Jacobian_Gauss_Legendre(Node_Coordinates)
    b_matrix = []
    for i in range(0, 4):
        result = np.dot(inv_jac[integration_point - 1], dh_drs[integration_point - 1][i])
        dh_dx, dh_dy = result[0], result[1]
        b_matrix.append([dh_dx, 0, 0, dh_dy, dh_dy, dh_dx])

    return b_matrix

print(B_Matrix_Assembly([[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]], 3)) # B matrix of node 1
'''