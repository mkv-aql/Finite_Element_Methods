__author__ = 'mkv-aql'
from sympy import Symbol, cos, factor, symbols, simplify
import math
import numpy as np

#Node coordinates = [Node 1, Node 2, Node 3, Node 4]
NC = [[0, 0], [1, 0], [1, 1], [0, 1]]
'''
##Jacobian calculation with Numerical Integration
#Linear shape function for 4 node element in 2D (N1 = h1)
r, s = symbols('r,s')
h1 = (1/4)*(1-r)*(1-s)
h2 = (1/4)*(1+r)*(1-s)
h3 = (1/4)*(1+r)*(1+s)
h4 = (1/4)*(1-r)*(1+s)

#Node coordinates = [Node 1, Node 2, Node 3, Node 4]
NC = [[0, 0], [1, 0], [1, 1], [0, 1]]

x_rs = h1*NC[0][0] + h2*NC[1][0] + h3*NC[2][0] + h4*NC[3][0]
y_rs = h1*NC[0][1] + h2*NC[1][1] + h3*NC[2][1] + h4*NC[3][1]

print(x_rs)
print(y_rs)

print("dx = ", simplify(x_rs))
print("dy = ", simplify(y_rs))

#Jacobian, J = [[dx/dr, dy/dr], [dx/ds, dy/ds]]
J = [[x_rs.diff(r), y_rs.diff(r)], [x_rs.diff(s), y_rs.diff(s)]]
print("Jacobian = ", J)

#Determinant of Jacobian
J_det = J[0][0]*J[1][1] - J[0][1]*J[1][0]
print("Det J = ", J_det)


print("")
print("")
print("")


##Jacobian calculation with Gauss Points or Quadrature Points
#Node coordinates = [Node 1, Node 2, Node 3, Node 4]
NC_1 = [[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]]

#Bilinear shape function for 4 node element in 2D (N1 = h1), dh1/dr = h1_r, dh1/ds = h1_s, and so on.
h1_r = h1.diff(r)
h1_s = h1.diff(s)
h2_r = h2.diff(r)
h2_s = h2.diff(s)
h3_r = h3.diff(r)
h3_s = h3.diff(s)
h4_r = h4.diff(r)
h4_s = h4.diff(s)

#substitute r and s values
l_rg = math.sqrt(1/3) #location of rg and sg
weight = 1 #Wg
combination_of_r_s = [[-(l_rg), -(l_rg)], [(l_rg), -(l_rg)], [(l_rg), (l_rg)], [-(l_rg), (l_rg)]]
value_1, value_2 = math.sqrt(1/3), math.sqrt(1/3)
h1_r = h1_r.subs({r:value_1, s:value_2})
h1_s = h1_s.subs({r:value_1, s:value_2})
h2_r = h2_r.subs({r:value_1, s:value_2})
h2_s = h2_s.subs({r:value_1, s:value_2})
h3_r = h3_r.subs({r:value_1, s:value_2})
h3_s = h3_s.subs({r:value_1, s:value_2})
h4_r = h4_r.subs({r:value_1, s:value_2})
h4_s = h4_s.subs({r:value_1, s:value_2})

#Jacobian, J = [[dx/dr, dy/dr], [dx/ds, dy/ds]] = [[A, B], [C, D]]
A = h1_r*NC_1[0][0] + h2_r*NC_1[1][0] + h3_r*NC_1[2][0] + h4_r*NC_1[3][0]
C = h1_s*NC_1[0][0] + h2_s*NC_1[1][0] + h3_s*NC_1[2][0] + h4_s*NC_1[3][0]
B = h1_r*NC_1[0][1] + h2_r*NC_1[1][1] + h3_r*NC_1[2][1] + h4_r*NC_1[3][1]
D = h1_s*NC_1[0][1] + h2_s*NC_1[1][1] + h3_s*NC_1[2][1] + h4_s*NC_1[3][1]

J_1 = [[A, B], [C, D]]
print("Jacobian_1 = ", J_1)

#Determinant of Jacobian
J_det_1 = A*D - B*C
print("Det J_1 = ", J_det_1)
'''

def Jacobian_NI(coordinate):
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

    x_rs_dr = x_rs.diff(r)
    y_rs_dr = y_rs.diff(r)
    x_rs_ds = x_rs.diff(s)
    y_rs_ds = y_rs.diff(s)

    # Jacobian, J = [[dx/dr, dy/dr], [dx/ds, dy/ds]]
    J = [[x_rs_dr, y_rs_dr], [x_rs_ds, y_rs_ds]]
    print("Jacobian = ", J)

    # Determinant of Jacobian
    J_det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    print("Det J = ", J_det)

    return J_det

print(Jacobian_NI(NC))

'''
for r_s in combination_of_r_s:
    print(r_s[0]) #value of r
    print(r_s[1]) #value of s
'''

def Jacobian_Gauss_Legendre(coordinate):
    # Linear shape function for 4 node element in 2D (N1 = h1)
    r, s = symbols('r,s')
    h1 = (1 / 4) * (1 - r) * (1 - s)
    h2 = (1 / 4) * (1 + r) * (1 - s)
    h3 = (1 / 4) * (1 + r) * (1 + s)
    h4 = (1 / 4) * (1 - r) * (1 + s)

    # Node coordinates = [Node 1, Node 2, Node 3, Node 4]
    NC = coordinate

    #Parameters for Gauss Quadrature
    l_rg = math.sqrt(1 / 3)  # location of rg and sg
    weight = 1  # Wg
    c_r_s = [[-(l_rg), -(l_rg)], [(l_rg), -(l_rg)], [(l_rg), (l_rg)], [-(l_rg), (l_rg)]] # Combination of r and s

    for r_s in c_r_s:
        print("")
        print("")
        #Iterate through the combination of r and s
        print("combination of r and s = ", r_s)
        index = c_r_s.index(r_s) # index of r_s in c_r_s

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

        # Jacobian, J = [[dx/dr, dy/dr], [dx/ds, dy/ds]] = [[A, B], [C, D]]
        dx_dr = dh1_dr * NC[0][0] + dh2_dr * NC[1][0] + dh3_dr * NC[2][0] + dh4_dr * NC[3][0]
        dx_ds = dh1_ds * NC[0][0] + dh2_ds * NC[1][0] + dh3_ds * NC[2][0] + dh4_ds * NC[3][0]
        dy_dr = dh1_dr * NC[0][1] + dh2_dr * NC[1][1] + dh3_dr * NC[2][1] + dh4_dr * NC[3][1]
        dy_ds = dh1_ds * NC[0][1] + dh2_ds * NC[1][1] + dh3_ds * NC[2][1] + dh4_ds * NC[3][1]

        Jac = [[dx_dr, dy_dr], [dx_ds, dy_ds]]
        print("Jacobian = ", Jac)

        # Determinant of Jacobian
        Jac_det = dx_dr * dy_ds - dx_ds * dy_dr
        print("Det J = ", Jac_det)

        # Convert to matrix
        Jac = np.array(Jac, dtype=float)

        # Inverse of Jacobian
        Jac_inv = np.linalg.inv(Jac)
        print("Jacobian Inverse = ", Jac_inv)



    return 0


print("")
print("")
Jacobian_Gauss_Legendre([[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]])
