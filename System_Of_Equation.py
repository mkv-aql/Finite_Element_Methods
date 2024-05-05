__author__ = 'mkv-aql'

import numpy as np

# 3x0 - x1 + 4x2 = 2,
# 17x0 + 2x1 + x2 = 14,
# x0 + 12x1 - 77x2 = 54
a = np.array([[3, -1, 4], [17, 2, 1], [1, 12, -77]])
b = np.array([2, 14, 54])

x = np.linalg.solve(a, b)
print(x)
np.allclose(np.dot(a, x), b) #Check if the solution is correct

c = np.array([[32.5, -1.16, -10], [-1.16, 12.5, 0], [-10, 0, 10]])
d = np.array([0, 0, 10])
y = np.linalg.solve(c, d)
print(y)

#Dataframes
import pandas as pd

df = pd.DataFrame(a, columns=['x0', 'x1', 'x2'])
df2 = pd.DataFrame(b)

df = df.to_numpy()
df2 = df2.to_numpy()

print(np.linalg.solve(df, df2))