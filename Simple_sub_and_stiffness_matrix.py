__author__ = 'mkv-aql'
import math

from sympy import Symbol, cos, factor, symbols

import numpy as np

# list of iteration combiantions
iteration_combinations = [[-(math.sqrt(1 / 3)), -(math.sqrt(1 / 3))],
                          [math.sqrt(1 / 3), -(math.sqrt(1 / 3))],
                          [math.sqrt(1 / 3), math.sqrt(1 / 3)],
                          [-(math.sqrt(1 / 3)), math.sqrt(1 / 3)]]

def calc_iteration(r_cor, s_cor):
    if (r_cor, s_cor) == (-(math.sqrt(1 / 3)), -(math.sqrt(1 / 3))):
        iteration = 1
    elif (r_cor, s_cor) == (math.sqrt(1 / 3), -(math.sqrt(1 / 3))):
        iteration = 2
    elif (r_cor, s_cor) == (math.sqrt(1 / 3), math.sqrt(1 / 3)):
        iteration = 3
    elif (r_cor, s_cor) == (-(math.sqrt(1 / 3)), math.sqrt(1 / 3)):
        iteration = 4
    else:
        iteration = None

    # weights
    wg = 1

    # Coordinates of the nodes
    natural_cor = [[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]]  # Parameter

    # Linear shape functions for quadrilateral element
    r, s = symbols('r,s')
    h1 = (1 / 4) * (1 - r) * (1 - s)
    h2 = (1 / 4) * (1 + r) * (1 - s)
    h3 = (1 / 4) * (1 + r) * (1 + s)
    h4 = (1 / 4) * (1 - r) * (1 + s)

    # Bilinear shape function for 4 node element in 2D (N1 = h1), dh1/dr = h1_r, dh1/ds = h1_s, and so on.
    dh1_dr = h1.diff(r)
    dh1_ds = h1.diff(s)
    dh2_dr = h2.diff(r)
    dh2_ds = h2.diff(s)
    dh3_dr = h3.diff(r)
    dh3_ds = h3.diff(s)
    dh4_dr = h4.diff(r)
    dh4_ds = h4.diff(s)

    # Insert r_cor and s_cor into the bilinear shape function
    dh1_dr = dh1_dr.subs({r: r_cor, s: s_cor})
    dh1_ds = dh1_ds.subs({r: r_cor, s: s_cor})
    dh2_dr = dh2_dr.subs({r: r_cor, s: s_cor})
    dh2_ds = dh2_ds.subs({r: r_cor, s: s_cor})
    dh3_dr = dh3_dr.subs({r: r_cor, s: s_cor})
    dh3_ds = dh3_ds.subs({r: r_cor, s: s_cor})
    dh4_dr = dh4_dr.subs({r: r_cor, s: s_cor})
    dh4_ds = dh4_ds.subs({r: r_cor, s: s_cor})

    # Insert into natural coordinates
    dx_dr = dh1_dr * natural_cor[0][0] + dh2_dr * natural_cor[1][0] + dh3_dr * natural_cor[2][0] + dh4_dr * \
            natural_cor[3][0]
    dx_ds = dh1_ds * natural_cor[0][0] + dh2_ds * natural_cor[1][0] + dh3_ds * natural_cor[2][0] + dh4_ds * \
            natural_cor[3][0]
    dy_dr = dh1_dr * natural_cor[0][1] + dh2_dr * natural_cor[1][1] + dh3_dr * natural_cor[2][1] + dh4_dr * \
            natural_cor[3][1]
    dy_ds = dh1_ds * natural_cor[0][1] + dh2_ds * natural_cor[1][1] + dh3_ds * natural_cor[2][1] + dh4_ds * \
            natural_cor[3][1]

    # Jacobian
    Jac = [[dx_dr, dy_dr], [dx_ds, dy_ds]]

    # Determinant of Jacobian
    det_Jac = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0]

    # Inverse matrix of Jacobian
    Jac_inv = np.array([[Jac[1][1], -Jac[0][1]], [-Jac[1][0], Jac[0][0]]])
    Jac_inv = Jac_inv * (1 / det_Jac)

    # Derivative of shape function for (x,y) or (x,y,z) if 3D
    [dh1_dx, dh1_dy] = np.dot(Jac_inv, [dh1_dr, dh1_ds])
    [dh2_dx, dh2_dy] = np.dot(Jac_inv, [dh2_dr, dh2_ds])
    [dh3_dx, dh3_dy] = np.dot(Jac_inv, [dh3_dr, dh3_ds])
    [dh4_dx, dh4_dy] = np.dot(Jac_inv, [dh4_dr, dh4_ds])

    # B matrix
    B1 = np.array([[dh1_dx, 0.0], [0.0, dh1_dy], [dh1_dy, dh1_dx]])
    B2 = np.array([[dh2_dx, 0.0], [0.0, dh2_dy], [dh2_dy, dh2_dx]])
    B3 = np.array([[dh3_dx, 0.0], [0.0, dh3_dy], [dh3_dy, dh3_dx]])
    B4 = np.array([[dh4_dx, 0.0], [0.0, dh4_dy], [dh4_dy, dh4_dx]])

    # Material matrix
    E = 10  # Small number just for simplifcation
    v = 0.3
    C = (E / (1 - v ** 2)) * np.array([[1, v, 0], [v, 1, 0], [0, 0, (1 - v) / 2]])

    # Sub Stiffness matrix
    B1_transpose = B1.transpose()
    # print("\n\nB1 transpose: \n", B1_transpose)
    B2_transpose = B2.transpose()
    B3_transpose = B3.transpose()
    B4_transpose = B4.transpose()

    sub_stiffness_1_1 = B1_transpose.dot(C).dot(B1) * det_Jac * wg
    sub_stiffness_2_1 = B1_transpose.dot(C).dot(B2) * det_Jac * wg
    sub_stiffness_3_1 = B1_transpose.dot(C).dot(B3) * det_Jac * wg
    sub_stiffness_4_1 = B1_transpose.dot(C).dot(B4) * det_Jac * wg
    sub_stiffness_2_2 = B2_transpose.dot(C).dot(B2) * det_Jac * wg
    sub_stiffness_3_2 = B2_transpose.dot(C).dot(B3) * det_Jac * wg
    sub_stiffness_4_2 = B2_transpose.dot(C).dot(B4) * det_Jac * wg
    sub_stiffness_3_3 = B3_transpose.dot(C).dot(B3) * det_Jac * wg
    sub_stiffness_4_3 = B3_transpose.dot(C).dot(B4) * det_Jac * wg
    sub_stiffness_4_4 = B4_transpose.dot(C).dot(B4) * det_Jac * wg

    # convert r_cor and s_cor to 4 decimal places
    # r_cor = round(r_cor, 3)
    # s_cor = round(s_cor, 3)

    print(f"\n\nSub Stiffness Matrix K(I,J) = (J=1,I=1) \nIteration: {iteration} = ({round(r_cor, 3)},{round(s_cor, 3)}) \n{sub_stiffness_1_1}")
    return sub_stiffness_1_1, sub_stiffness_2_1, sub_stiffness_3_1, sub_stiffness_4_1, sub_stiffness_2_2, sub_stiffness_3_2, sub_stiffness_4_2, sub_stiffness_3_3, sub_stiffness_4_3, sub_stiffness_4_4

sum = np.zeros((2, 2))
# loop through the iteration_combinations
for i in iteration_combinations:
    # print("i = ", i)
    r_cor = i[0]
    s_cor = i[1]

    # sub_stiffness_1_1, _, _, _ = calc_iteration(r_cor, s_cor)
    # _, sub_stiffness_2_2, _, _ = calc_iteration(r_cor, s_cor)
    # _, _, sub_stiffness_3_3, _ = calc_iteration(r_cor, s_cor)
    # _, _, _, sub_stiffness_4_4 = calc_iteration(r_cor, s_cor)
    sub_stiffness_1_1, sub_stiffness_2_1, sub_stiffness_3_1, sub_stiffness_4_1, sub_stiffness_2_2, sub_stiffness_3_2, sub_stiffness_4_2, sub_stiffness_3_3, sub_stiffness_4_3, sub_stiffness_4_4 = calc_iteration(r_cor, s_cor)
    # sum = np.add(sum, calc_iteration(r_cor, s_cor))
    sum = np.add(sum, sub_stiffness_4_3) # to get K43 of all iteration_combinations


print(f'\n\nSub Stiffness matrix (J,I) = (4,3) \nAll Iteration (summed up) \n{sum}')

