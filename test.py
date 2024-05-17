import numpy as np


def calculate_nodes(start, end, h):
    # calculate points :
    nodes = np.arange(start, end + h - 0.000001, h)

    return nodes


def u1(x, y):
    return 4 * x ** 4 + 8 * y ** 4


def f1(x, y):
    return 48 * (x ** 2 + 2 * y ** 2)


x_0 = 0
x_n = 0.5
y_0 = 0
y_n = 1


h_x = 0.1
h_y = 0.1

nx = int((x_n - x_0) / h_x + 1)
ny = int((y_n - y_0) / h_y + 1)

x_nodes = calculate_nodes(x_0, x_n, h_x)
y_nodes = calculate_nodes(y_0, y_n, h_y)

exact_solution = np.zeros((nx, ny))

for i in range(nx):
    for j in range(ny):
        exact_solution[i, j] = u1(x_nodes[i], y_nodes[j])

boundary = np.zeros((nx, ny))
boundary[0, :] = exact_solution[0, :]
boundary[-1, :] = exact_solution[-1, :]
boundary[:, 0] = exact_solution[:, 0]
boundary[:, -1] = exact_solution[:, -1]

rhs = np.zeros((nx, ny))

for i in range(1, nx-1):
    for j in range(1, ny-1):
        rhs[i, j] = f1(x_nodes[i], y_nodes[j])

approx = np.zeros((nx, ny))

U = boundary

for i in range(1, nx-1):
    for j in range(1, ny-1):
        approx[i, j] = (U[i-1, j] + U[nx-1, j] + U[i, j-1] + U[i, ny-1] - h_x ** 2 * rhs[i, j])/4

U += approx
for t in range(50):
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            approx[i, j] = (U[i-1, j] + U[i+1, j] + U[i, j-1] + U[i, j+1] - h_x ** 2 * rhs[i, j])/4
            U[i, j] = approx[i, j]


print(np.abs(exact_solution-U).max())




