__author__ = 'mkv-aql'

from sympy import Symbol, cos, factor, symbols, integrate

import numpy as np

# Linear shape functions for quadrilateral element, h = high, l = low
x, y, x_l, x_h, y_l, y_h, inv_area, v, E = symbols('x, y, x_l, x_h, y_l, y_h, 1/A, v, E')
h1 = (inv_area)*(x_h - x)*(y_h - y)
h2 = (inv_area)*(x - x_l)*(y_h - y)
h3 = (inv_area)*(x - x_l)*(y - y_l)
h4 = (inv_area)*(x_h - x)*(y - y_l)



# Parameters
to = -500 # Parameter edge load
area = 15   # Parameter, example 1
# area = 1 # Parameter, example 2
results = h2 * to
cor = [[1, 1], [6, 1], [6, 4], [1, 4]] # Parameter, example 1
# cor = [[0, 0], [1.5, 0], [1.5, 2/3], [0, 2/3]] # Parameter, example 2

edge_load_cor = [[6, 1], [6, 4]] # Parameter, along y axis, example 1
# edge_load_cor = [[3/2, 0], [3/2, 2/3]] # Parameter, along y-axis, example 2

print(f'\nResults: {results}')

results = results.subs({inv_area: 1/area, x_l: cor[0][0], x_h: cor[1][0], y_l: cor[0][1], y_h: cor[2][1] , x: edge_load_cor[0][0]}) # for edge load towards y-axis
# results = results.subs({inv_area: 1/area, x_l: cor[0][0], x_h: cor[1][0], y_l: cor[0][1], y_h: cor[2][1] , y: edge_load_cor[1][1]}) # for edge load towards x-axis

print(f'\nResults: {results}')


# results = integrate(results, (x, cor[0][0], cor[1][0]), (y, cor[0][1], cor[2][1]))
# results = integrate(results, (x, cor[0][0], cor[1][0])) # integrate only x
results = integrate(results, (y, cor[0][1], cor[2][1])) # integrate only y

print(f'\nResults: {results}')


