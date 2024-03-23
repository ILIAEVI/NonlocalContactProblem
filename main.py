"""
d^2u_1/dx^2=f_1(x)  where 0 < x < c
d^2u_2/dx^2=f_2(x)  where c < x < 1

first boundary condition for every iteration:
u_1(0) = 0
u_2(1) = 0

Non local boundary condition:  iteration is k,
we denote iteration as degree ^(k)
U1[c]^(k) = U2^[c](k) = a * U1(x_1)^(k-1) + b * U2(x_2)^(k-1) + f_0

where a = const > 0  and b = const > 0, a + b < 1,
and we take also f_0

with this iteration process we reduce problem to the solution of the
sequence of classical dirichlet problem


example:
d^2u_1/dx^2=2x^2  where 0 < x < 1/2
d^2u_2/dx^2=2x+5  where 1/2 < x < 1

u_1(1/2)=u_2(1/2)=2

"""

import numpy as np


def compute_nodal_points(h, x_0, x_n):
    nodal_points = np.arange(x_0, x_n + h, h)
    return nodal_points


def compute_matrix_for_difference_method(size):
    matrix = np.zeros((size, size))
    matrix[0, 0] = 1
    matrix[size-1, size-1] = 1
    for i in range(1, size-1):
        matrix[i, i - 1] = 1
        matrix[i, i] = -2
        matrix[i, i + 1] = 1
    return matrix


def left_finite_difference_method(g, nodal_points, f, h):
    n = len(nodal_points)
    y = np.zeros(n)

    # boundary condition u(0)=0 and u(0.5)=u(1/2)
    y[0] = 0
    y[n-1] = g[2]

    matrix = compute_matrix_for_difference_method(n)
    for i in range(1, n-1):
        y[i] = f(nodal_points[i]) * h ** 2

    result = np.linalg.solve(matrix, y)

    """
    F(x_i) = (y(x_i+1) - 2 * y(x_i) + y(x_i-1))/ h^2
    """
    return result


def right_finite_difference_method(g, nodal_points, f, h):
    n = len(nodal_points)
    y = np.zeros(n)

    # boundary condition u(0.5)=u[2] and u(1)=0
    y[0] = g[2]
    y[n-1] = 0

    matrix = compute_matrix_for_difference_method(n)
    for i in range(1, n - 1):
        y[i] = f(nodal_points[i]) * h ** 2

    result = np.linalg.solve(matrix, y)

    """
    F(x_i) = (y(x_i+1) - 2 * y(x_i) + y(x_i-1))/ h^2
    """

    return result


def non_local_contact_problem_iteration(nodes, u_1, u_2, f_1, f_2, k, h):
    # define gama_1 and gama_2 for nonlocal contact boundary condition
    # gama_1 + gama_2 < 1, gama_1 > 0 and gama_2 > 0
    gama_1 = 0.25
    gama_2 = 0.25
    # n must be 5, nodes= [0, 0.25, 0.5, 0.75, 1]
    n = len(nodes)
    # f_0 is necessary to make sure that the contact condition is correct
    f_0 = u_1(nodes[2]) - (gama_1 * u_1(nodes[1]) + gama_2 * u_2(nodes[3]))
    left_nodes = compute_nodal_points(h, nodes[0], nodes[2])
    right_nodes = compute_nodal_points(h, nodes[2], nodes[4])
    if n != 5:
        raise print('please check main nodes, length of nodes is not 5.'
                    'make sure that nodes length is five or change iteration code')

    u = np.zeros(n)
    left_mid_point = int((len(left_nodes)-1)/2)
    right_mid_point = int((len(right_nodes)-1)/2)
    for i in range(k):
        u[2] = gama_1 * u[1] + gama_2 * u[3] + f_0
        y_left = left_finite_difference_method(u, left_nodes, f_1, h)
        y_right = right_finite_difference_method(u, right_nodes, f_2, h)
        u[1] = y_left[left_mid_point]
        u[3] = y_right[right_mid_point]
        print(f"Iteration {i+1}: u(1/2) = {u[2]}")
        print(f"error: {abs(u_1(nodes[2]) - u[2])}")
        print('-'*50)
        print(f"iteration: {i+1}, u(1/4): {u[1]}")
        print(f"iteration:{i+1}, error u(1/4): {abs(u_1(nodes[1])-u[1])}")
        print('-' * 50)
        print(f"iteration: {i+1}, u(3/4): {u[3]}")
        print(f"iteration:{i + 1}, error u(3/4): {abs(u_2(nodes[3]) - u[3])}")
        print('*' * 50)


def calculate_absolute_error(u_1, points, k, apr):
    exact_result = u_1(points[2])  # u(1/2)
    errors = {}
    for i in range(k):
        errors[i+1] = abs(exact_result - apr[2])
        print(f"Absolute Error on {i+1} iteration: {errors[i+1]}")
    minimum_error = min(errors, key=errors.get)
    print(f"Iteration: {minimum_error} We get minimum error: {min(errors.values())}")
    print(f"fExact value in 1/2: {exact_result}")


if __name__ == "__main__":
    def u1(x):
        return 1/48*(8 * x**4 + 71 * x)

    def u2(x):
        return 1 / 6 * (2 * x ** 3 + 15 * x ** 2 - 35 * x + 18)

    def f1(x):
        return 2 * x**2

    def f2(x):
        return 2 * x + 5

    initial_nodal_points = np.arange(0, 1.25, 0.25)
    # iteration number:
    iteration = 14
    # define division number
    division = 0.001

    non_local_contact_problem_iteration(
        initial_nodal_points,
        u1, u2, f1, f2,
        iteration,
        division
    )

    # calculate_absolute_error(u1, initial_nodal_points, iteration_number, result)
