__author__ = 'mkv-aql'
import numpy as np
from fractions import Fraction
import sympy as sp
import pandas as pd

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