__author__ = 'mkv-aql'
import pandas as pd
import numpy as np

# Sample matrices with variables a, b, c, and d
data1 = np.array([[1, 2], [3, 4]])  # Variables a, b
data2 = np.array([[5, 6], [7, 8]])  # Variables b, c
data3 = np.array([[9, 10], [11, 12]])  # Variables c, d

# Convert to pandas DataFrames with column labels
df1 = pd.DataFrame(data1, columns=['a', 'b'])
df2 = pd.DataFrame(data2, columns=['b', 'c'])
df3 = pd.DataFrame(data3, columns=['c', 'd'])

# Initialize an empty DataFrame to combine everything
combined_df = pd.DataFrame()

# List of all DataFrames for iteration
dfs = [df1, df2, df3]

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

print(combined_df)