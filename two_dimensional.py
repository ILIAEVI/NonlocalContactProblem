import numpy as np


def u1(x1, y1):
    return 4*x1**4 + 8*y1**4


def u2(x1, y1):
    return 2 * x1 ** 3 + 8 * y1 ** 4


def f1(x1, y1):
    return 48 * (x1 ** 2 + 2 * y1 ** 2)


def f2(x1, y1):
    return 12 * (x1 + 8 * y1 ** 2)

# ____________________________________________


gama_1 = 0.25
gama_2 = 0.25

x_main_points = np.arange(0, 1.25, 0.25)
n = len(x_main_points)


if n != 5:
    raise print('please check main nodes, length of nodes is not 5.'
                'make sure that nodes length is five or change iteration code')

u = np.zeros((n, n))


def f_0(y1):
    return u1(x_main_points[2], y1) - (gama_1 * u1(x_main_points[1], y1) + gama_2 * u2(x_main_points[3], y1))


# y is from  0 to 1
y_0 = 0
y_n = 1
# note that if u take division number more than 0.1, it may result in significant errors.
# request: division number <= 0.1
h_x = 0.1  # x axis division
h_y = 0.1  # y axis division



left_x_nodes = np.arange(0, 1 + h_x-0.000001, h_x)
left_y_nodes = np.arange(y_0, y_n + h_y-0.000001, h_y)
left_X, left_Y = np.meshgrid(left_x_nodes, left_y_nodes)

right_x_nodes = np.arange(x_main_points[2], x_main_points[4] + h_x-0.000001, h_x)
right_y_nodes = np.arange(y_0, y_n + h_y, h_y)

right_X, right_Y = np.meshgrid(right_x_nodes, right_y_nodes)

left_x_mid_point = int((len(left_x_nodes)-1)/2)
right_x_mid_point = int((len(right_x_nodes)-1)/2)

# calculate exact solution for every point
exact_u_left = np.zeros((len(left_x_nodes), len(left_y_nodes)))
exact_u_right = np.zeros((len(right_x_nodes), len(right_y_nodes)))


for i, x in enumerate(left_x_nodes):
    for j, y in enumerate(left_y_nodes):
        exact_u_left[i, j] = u1(x, y)

for i, x in enumerate(right_x_nodes):
    for j, y in enumerate(right_y_nodes):
        exact_u_right[i, j] = u2(x, y)

full_exact_u = np.concatenate((exact_u_left[:-1, :], exact_u_right), axis=0)


# build matrix
n_x = len(left_x_nodes)
n_y = len(left_y_nodes)
print(n_x)
N = n_x * n_y
matrix = np.zeros((N, N))

# diagonal

for i in range(0, n_x):
    for j in range(0, n_y):
        k = i + n_x*j
        matrix[k, k] = -4

# Lower Diagonal
for i in range(1, n_x):
    for j in range(0, n_y):
        matrix[i+n_x*j, i+n_x*j-1] = 1
# Upper Diagonal
for i in range(0, n_x-1):
    for j in range(0, n_y):
        matrix[i+n_x*j, i+n_x*j+1] = 1

# Lower I Matrix
for i in range(0, n_x):
    for j in range(1, n_y):
        matrix[i+n_x*j, i+n_x*(j-1)] = 1

# Upper I Matrix:
for i in range(0, n_x):
    for j in range(0, n_y-1):
        matrix[i+n_x*j, i+n_x*(j+1)] = 1

inv = np.linalg.inv(matrix)

# Right-hand side function for left problem
def rhs_left(t_1, t_2):
    return 48 * (t_1**2 + 2 * t_2**2)

# calculate right hand side values
# Calculate right-hand side values
rhs_left_values = np.zeros(N)

for i in range(n_x):
    for j in range(n_y):
        rhs_left_values[i + n_x * j] = f1(left_x_nodes[i + 1], left_y_nodes[j + 1]) * h_x**2


boundary = np.zeros((n_x+1, n_y+1))
boundary[0, :] = u1(left_x_nodes, 0)  # Bottom boundary
boundary[-1, :] = u1(left_x_nodes, 1)  # Top boundary
boundary[:, 0] = u1(0, left_y_nodes)  # Left boundary
boundary[:, -1] = u1(1, left_y_nodes)  # Right boundary



boundary = boundary.flatten()
print('bounaru new: \n', boundary)
print(len(boundary))
'''
rhs_left_values[0] = exact_u_left[0]
rhs_left_values[-1] = exact_u_left[-1]
rhs_left_values[:, 0] = exact_u_left[:, 0]
rhs_left_values[:, -1] = exact_u_left[:, -1]'''


#rhs_left_values = np.reshape(rhs_left_values, (-1, 1))
rhs_left_values = rhs_left_values.flatten()
print('rhs: \n', rhs_left_values)
#print('matrix: \n', matrix)

solution = np.linalg.solve(matrix, rhs_left_values-boundary)
solution = np.dot(inv, rhs_left_values-boundary)
print('solution: \n', solution)
#exact_u_left = np.reshape(exact_u_left, (-1, 1))
exact_u_left = exact_u_left.flatten()

print('exact: \n', exact_u_left)

