__author__ = 'mkv-aql'

import numpy as np
import math

gauss_quadrature = {1: {'weight': 2, 'location': [0]},
                    2: {'weight': 1, 'location': [-1/np.sqrt(3), 1/np.sqrt(3)]},
                    3: {'weight': 5/9, 'location': [-np.sqrt(3/5), 0, np.sqrt(3/5)]}
                    }

gauss_quadrature_r_s_combination = {1: [[0, 0]],
                                    2: [[-1/np.sqrt(3), -1/np.sqrt(3)],
                                        [1/np.sqrt(3), -1/np.sqrt(3)],
                                        [1/np.sqrt(3), 1/np.sqrt(3)],
                                        [-1/np.sqrt(3), 1/np.sqrt(3)]]
                                    }

'''
print(gauss_quadrature_r_s_combination[2])
for i in gauss_quadrature_r_s_combination[2]:
    print(i[0])
    print(i[1])
    print("")
    print("")
'''
a = gauss_quadrature_r_s_combination[2]
print(a)
l_rg = math.sqrt(1 / 3)  # location of rg and sg
weight = 1  # Wg
c_r_s = [[-(l_rg), -(l_rg)], [(l_rg), -(l_rg)], [(l_rg), (l_rg)], [-(l_rg), (l_rg)]]  # Combination of r and s
print(c_r_s)
