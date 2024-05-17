import numpy as np


def calculate_nodes(start, end, h):
    # calculate points :
    nodes = np.arange(start, end + h - 0.000001, h)

    return nodes


def calculate_exact_solution(nx, ny, u, x_points, y_points):
    # exact solution :
    ex_solution = np.zeros((nx, ny))

    for i in range(nx):
        for j in range(ny):
            ex_solution[i, j] = u(x_points[i], y_points[j])

    return ex_solution


# boundary:


def set_boundary_conditions(nx, ny, exact_sol):
    boundary = np.zeros((nx, ny))
    boundary[0, :] = exact_sol[0, :]
    boundary[-1, :] = exact_sol[-1, :]
    boundary[:, 0] = exact_sol[:, 0]
    boundary[:, -1] = exact_sol[:, -1]

    return boundary


def build_rhs_vector(nx, ny, x_points, y_points, boundary_conditions, h_1, h_2, f):

    rhs = np.zeros((nx, ny))

    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            rhs[i, j] = f(x_points[i], y_points[j]) * h_1 * h_2

    rhs += boundary_conditions

    return rhs.flatten()


def build_matrix_for_finite_difference_2d(nx, ny):
    n = nx * ny
    matrix = np.zeros((n, n))

    for i in range(0, nx + 1):
        matrix[i, i] = 1

    for j in range(2, nx):
        for i in range((j - 1) * nx + 1, j * nx - 1):
            matrix[i, i] = -4
            matrix[i, i + 1] = 1
            matrix[i, i - 1] = 1
            matrix[i, i - nx] = 1
            matrix[i, i + nx] = 1

        for i in range(j * nx - 1, j * nx + 1):
            matrix[i, i] = 1

    for i in range(n - nx, n):
        matrix[i, i] = 1

    return matrix


def max_absolute_error(approx, ex_sol):
    return max(abs(ex_sol.flatten()-approx))


if __name__ == "__main__":
    def u1(x, y):
        return 4 * x ** 4 + 8 * y ** 4


    def f1(x, y):
        return 48 * (x ** 2 + 2 * y ** 2)


    def u2(x, y):
        return 2 * x ** 3 + 8 * y ** 4


    def f2(x, y):
        return 12 * (x + 8 * y ** 2)

    # define domain
    x = np.arange(0, 1.25-0.00001, 0.25)
    x_0 = 0
    x_n = 1
    y_0 = 0
    y_n = 1

    # define division:
    h_x = 0.01
    h_y = 0.01

    # calculate grid number
    x_grid = int((x_n-x_0) / h_x + 1)
    y_grid = int((y_n-y_0) / h_y + 1)

    x_nodes = calculate_nodes(x_0, x_n, h_x)
    y_nodes = calculate_nodes(y_0, y_n, h_y)

    exact_solution = calculate_exact_solution(x_grid, y_grid, u1, x_nodes, y_nodes)

    bound = set_boundary_conditions(x_grid, y_grid, exact_solution)

    r = build_rhs_vector(x_grid, y_grid, x_nodes, y_nodes, bound, h_x, h_y, f1)

    A = build_matrix_for_finite_difference_2d(x_grid, y_grid)

    # print(A)
    approximation = np.linalg.solve(A, r)
    print(approximation)
    error = max_absolute_error(approximation, exact_solution)
    # print(exact_solution)
    print(error)
