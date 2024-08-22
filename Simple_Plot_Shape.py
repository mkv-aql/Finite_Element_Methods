__author__ = 'mkv-aql'
import matplotlib.pyplot as plt

# Step 1: Define the coordinates of the vertices of the shape
# In this example, let's create a simple quadrilateral
# cor = [[1, 1], [6, 1], [6, 4], [1, 4]] # Parameter, example 1
# cor = [[0, 0], [1.5, 0], [1.5, 2/3], [0, 2/3]] # Parameter, example 2
cor = [[-1, -1], [0.8, -1.2], [0.6, 1], [-1.2, 0.6]]  # Parameter, example 3
x = [1, 4, 4, 1, 1]  # x coordinates of the vertices, repeated the first vertex to close the shape
y = [1, 1, 3, 3, 1]  # y coordinates of the vertices
x = [cor[0][0], cor[1][0], cor[2][0], cor[3][0], cor[0][0]]  # x coordinates of the vertices
y = [cor[0][1], cor[1][1], cor[2][1], cor[3][1], cor[0][1]]  # y coordinates of the vertices

# Step 2: Create the plot
plt.plot(x, y, 'bo-', label='Shape')  # 'bo-' indicates blue color, circle markers, and solid line

# Step 3: Set up labels and title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Simple Quadrilateral')
plt.legend()

# Step 4: Display the plot
plt.grid(True)
plt.axis('equal')  # Ensure the aspect ratio is equal so the shape isn't distorted
plt.show()