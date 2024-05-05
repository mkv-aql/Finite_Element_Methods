__author__ = 'mkv-aql'
boundaries_dataframe_1 = ['u1', 'v1', 'u3', 'v3', 'v4']
Initial_DOF = ['u1', 'v1', 'u2', 'v2', 'u3', 'v3', 'u4', 'v4'] #Initial DOF of system
unique = set(Initial_DOF) ^ set(boundaries_dataframe_1) #Remove boundary conditions from DOF
Final_DOF = [x for x in Initial_DOF if x in unique] #Final DOF after removing boundary conditions
print(Final_DOF)