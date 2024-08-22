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


#Symbols
from sympy import symbols, Eq, solve

x, y, z = symbols('x y z')
eq1 = Eq(2*x + 3*y - z, 1)
eq2 = Eq(-x + 7*y + 2*z, -3)
eq3 = Eq(4*x - 2*y + 5*z, 2)

solution = solve((eq1, eq2, eq3), (x, y, z))
print(solution)

'''
If 1 or more values in x array is known: eg. x2 = 5
x = [x0, x1, 5]

System must be reduced ba removing x2 row and column

'''

import sympy as sp

# Define
# the
# coefficients
# of
# the
# equations
# Equation 1: 2x + 3y - z = 1
# Equation 2: -x + 7y + 2z = -3
# Equation 3: 4x - 2y + 5z = 2

coefficients = [
    [2, 3, -1, 1],  # coefficients of x, y, z, and the constant term for eq1
    [-1, 7, 2, -3],  # coefficients of x, y, z, and the constant term for eq2
    [4, -2, 5, 2]  # coefficients of x, y, z, and the constant term for eq3
]

# Known value for z
known_values = {'z': -2}


# Function to substitute known values into equations
def substitute_and_reduce(coefficients, known_values):
    reduced_coefficients = []

    for eq in coefficients:
        reduced_eq = eq[:-1]  # Copy the coefficients of x, y, z (without the constant term)
        constant_term = eq[-1]  # Copy the constant term

        # Substitute the known values
        for i, var in enumerate(['x', 'y', 'z']):
            if var in known_values:
                constant_term -= eq[i] * known_values[var]
                reduced_eq[i] = 0  # Remove the variable's contribution

        # Append the new reduced equation with the updated constant term
        reduced_eq.append(constant_term)
        reduced_coefficients.append(reduced_eq)

    return reduced_coefficients


# Reduce the system of equations
reduced_coefficients = substitute_and_reduce(coefficients, known_values)


# Function to solve a reduced 2x2 system of linear equations manually
def solve_2x2(coefficients):
    a1, b1, c1 = coefficients[0][:3]
    a2, b2, c2 = coefficients[1][:3]

    determinant = a1 * b2 - a2 * b1
    if determinant == 0:
        raise ValueError("The system has no unique solution.")

    x = (c1 * b2 - c2 * b1) / determinant
    y = (a1 * c2 - a2 * c1) / determinant

    return {'x': x, 'y': y}


# Remove the row where z is entirely zero (no variable contribution)
reduced_coefficients = [eq for eq in reduced_coefficients if any(eq[:-1])]

# Solve the remaining system of equations
solution = solve_2x2(reduced_coefficients)

# Add the known value back to the solution
solution.update(known_values)

# Print the reduced system of equations
print("Reduced system of equations:")
for eq in reduced_coefficients:
    print(f"{eq[0]}*x + {eq[1]}*y = {eq[2]}")

# Print the solution
print("\nSolution:")
for var, val in solution.items():
    print(f"{var} = {val}")