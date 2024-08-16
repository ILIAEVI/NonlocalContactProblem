import time

import numpy as np
import matplotlib.pyplot as plt
import csv
from concurrent.futures import ProcessPoolExecutor

def compute_nodes(start, end, h):
    # calculate points :
    nodes = np.arange(start, end + h - 0.000001, h)

    return nodes


def compute_exact_solution_matrix(u, x_grid, y_grid, x_points, y_points):
    exact_solution = np.zeros((x_grid, y_grid))

    for k in range(x_grid):
        for m in range(y_grid):
            exact_solution[k, m] = u(x_points[k], y_points[m])
    return exact_solution


def compute_boundary_condition(exact_sol, x_gird, y_grid):
    # There we calculate boundary conditions from exact solution:
    # In this case, we know the exact solution.
    boundary = np.zeros((x_gird, y_grid))
    boundary[0, :] = exact_sol[0, :]
    boundary[-1, :] = exact_sol[-1, :]
    boundary[:, 0] = exact_sol[:, 0]
    boundary[:, -1] = exact_sol[:, -1]

    return boundary


def rhs_vector(x_grid, y_grid, x_points, y_points, f):
    rhs = np.zeros((x_grid, y_grid))

    for k in range(1, x_grid - 1):
        for m in range(1, y_grid - 1):
            rhs[k, m] = f(x_points[k], y_points[m])

    return rhs


def fdm_iteration_left(x_grid, y_grid, boundary, contact_boundary, rhs, x_h, y_h, fdm_iter):
    approx = np.zeros((x_grid, y_grid))
    initial = boundary
    initial[-1] = contact_boundary[2]
    # Compute in advance to speed up iteration
    tau = 1 / (2 * (1 / x_h ** 2 + 1 / y_h ** 2))
    hx_2 = 1 / x_h ** 2
    hy_2 = 1 / y_h ** 2

    for k in range(1, x_grid - 1):
        for m in range(1, y_grid - 1):
            approx[k, m] = tau * (hx_2 * (initial[x_grid - 1, m] + initial[k - 1, m]) + hy_2 * (
                        initial[k, m - 1] + initial[k, y_grid - 1]) - rhs[k, m])

    approx += initial

    for t in range(fdm_iter):
        for k in range(1, x_grid - 1):
            for m in range(1, y_grid - 1):
                approx[k, m] = tau * (hx_2 * (approx[k + 1, m] + approx[k - 1, m]) + hy_2 * (
                            approx[k, m + 1] + approx[k, m - 1]) - rhs[k, m])

    return approx


def fdm_iteration_right(x_grid, y_grid, boundary, contact_boundary, rhs, x_h, y_h, fdm_iter):
    approx = np.zeros((x_grid, y_grid))
    initial = boundary
    initial[0] = contact_boundary[2]
    # Compute in advance to speed up iteration
    tau = 1 / (2 * (1 / x_h ** 2 + 1 / y_h ** 2))
    hx_2 = 1 / x_h ** 2
    hy_2 = 1 / y_h ** 2

    for k in range(1, x_grid - 1):
        for m in range(1, y_grid - 1):
            approx[k, m] = tau * (hx_2 * (initial[x_grid - 1, m] + initial[k - 1, m]) + hy_2 * (
                        initial[k, m - 1] + initial[k, y_grid - 1]) - rhs[k, m])

    approx += initial

    for t in range(fdm_iter):
        for k in range(1, x_grid - 1):
            for m in range(1, y_grid - 1):
                approx[k, m] = tau * (hx_2 * (approx[k + 1, m] + approx[k - 1, m]) + hy_2 * (
                            approx[k, m + 1] + approx[k, m - 1]) - rhs[k, m])

    return approx


def update_plot_3d(left_x_nodes, right_x_nodes, y_nodes, exact_solution, approximation, k, xh, yh):
    x_nodes = np.concatenate((left_x_nodes[:-1], right_x_nodes))
    X, Y = np.meshgrid(x_nodes, y_nodes)

    fig = plt.figure(figsize=(14, 7))

    ax1 = fig.add_subplot(121, projection='3d')
    ax1.plot_surface(X, Y, exact_solution.T, cmap='viridis')
    ax1.set_title('Exact Solution')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('u')

    ax2 = fig.add_subplot(122, projection='3d')
    ax2.plot_surface(X, Y, approximation.T, cmap='inferno')
    ax2.set_title(f'Approximation at iteration {k + 1}')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('u')

    plt.tight_layout()
    plt.show()


def plot_surface(left_x_nodes, right_x_nodes, y_nodes, exact_solution, approximation, k, h_1, h_2):
    x_nodes = np.concatenate((left_x_nodes[:-1], right_x_nodes))
    X, Y = np.meshgrid(x_nodes, y_nodes)

    fig = plt.figure(figsize=(14, 7))

    ax1 = fig.add_subplot(121, projection='3d')
    surf1 = ax1.plot_surface(X, Y, exact_solution.T, cmap='viridis', edgecolor='k', linewidth=0.5)
    fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10, pad=0.1)
    ax1.set_title('Exact Solution')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('u')
    ax1.view_init(elev=25, azim=65)  # Adjust the elevation and azimuth for better view

    ax2 = fig.add_subplot(122, projection='3d')
    surf2 = ax2.plot_surface(X, Y, approximation.T, cmap='viridis', edgecolor='k', linewidth=0.5)
    fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10, pad=0.1)
    ax2.set_title(f'Approximation at iteration {k + 1}, \n $h_{{x}} = {h_1},\ h_{{y}} = {h_2}$')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('u')
    ax2.view_init(elev=25, azim=65)  # Adjust the elevation and azimuth for better view

    plt.tight_layout()
    plt.show()


def plot_error(left_x_nodes, right_x_nodes, y_nodes, errors, k, h_1, h_2):
    x_nodes = np.concatenate((left_x_nodes[:-1], right_x_nodes))
    X, Y = np.meshgrid(x_nodes, y_nodes)

    fig = plt.figure(figsize=(12, 6))

    ax2 = fig.add_subplot(111, projection='3d')
    surf2 = ax2.plot_surface(X, Y, errors.T, cmap='viridis', edgecolor='k', linewidth=0.5)
    fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10, pad=0.1)
    ax2.set_title(f'Approximation Error at iteration {k + 1}, \n $h_{{x}} = {h_1},\ h_{{y}} = {h_2}$')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.view_init(elev=25, azim=65)

    plt.tight_layout()
    plt.show()


def plot_iteration_errors(iteration_errors, h_1, h_2):
    iterations, errors = zip(*iteration_errors)
    plt.figure(figsize=(10, 6))
    plt.plot(iterations, errors, label='Approximation Error')
    plt.title(f'Non-local Contact Problem Approximation.    $h_{{x}} = {h_1},\ h_{{y}} = {h_2}$')
    plt.yscale('log')
    plt.xlabel('Iteration Number')
    plt.ylabel('Max Absolute Error')
    plt.legend()
    plt.grid(True)
    plt.axhline(0, color='black', linewidth=1)
    plt.axvline(0, color='black', linewidth=1)
    plt.show()


def non_local_contact_problem_iter(initial_x_points, initial_y_points, u_1, u_2, f_1, f_2, iteration, fdm_iter, xh,
                                   yh, ):
    # define gama_1 and gama_2 for nonlocal contact boundary condition
    # gama_1 + gama_2 <= 1, gama_1 > 0 and gama_2 > 0
    gama_1 = 0.5
    gama_2 = 0.25

    # n must be 5, nodes= [0, 0.25, 0.5, 0.75, 1]
    n_1 = len(initial_x_points)
    n_2 = len(initial_y_points)

    # grid
    left_nx = int((initial_x_points[2] - initial_x_points[0]) / xh + 1)
    right_nx = int((initial_x_points[4] - initial_x_points[2]) / xh + 1)
    nx = int((initial_x_points[4] - initial_x_points[0]) / xh + 1)
    ny = int((initial_y_points[1] - initial_y_points[0]) / yh + 1)

    if n_2 != 2:
        raise print("Check initial y points")

    if n_1 != 5:
        raise print('please check main nodes, length of nodes is not 5.'
                    'make sure that nodes length is five or change iteration code')

    left_x_nodes = compute_nodes(initial_x_points[0], initial_x_points[2], xh)
    right_x_nodes = compute_nodes(initial_x_points[2], initial_x_points[4], xh)
    y_nodes = compute_nodes(initial_y_points[0], initial_y_points[1], yh)

    left_x_mid_point = int((len(left_x_nodes) - 1) / 2)
    right_x_mid_point = int((len(right_x_nodes) - 1) / 2)

    # calculate exact solution in every node.
    u_left = compute_exact_solution_matrix(u_1, left_nx, ny, left_x_nodes, y_nodes)
    u_right = compute_exact_solution_matrix(u_2, right_nx, ny, right_x_nodes, y_nodes)

    exact_solution = np.concatenate((u_left[:-1, :], u_right), axis=0)

    # calculate f_0 vector for contact boundary condition
    f_0 = np.zeros(ny)
    for i in range(ny):
        f_0[i] = u_1(initial_x_points[2], y_nodes[i]) - (
                    gama_1 * u_1(initial_x_points[1], y_nodes[i]) + gama_2 * u_2(initial_x_points[3], y_nodes[i]))

    v = np.zeros((n_1, ny))
    # in my example i dont need compute boundary, while it is 0 in every node
    #left_boundaries = compute_boundary_condition(u_left, left_nx, ny)
    #right_boundaries = compute_boundary_condition(u_right, right_nx, ny)
    left_boundaries = np.zeros((left_nx, ny))
    right_boundaries = np.zeros((right_nx, ny))

    left_rhs = rhs_vector(left_nx, ny, left_x_nodes, y_nodes, f_1)
    right_rhs = rhs_vector(right_nx, ny, right_x_nodes, y_nodes, f_2)

    errors = []
    iteration_errors = []

    for i in range(iteration):
        v[2] = gama_1 * v[1] + gama_2 * v[3] + f_0
        with ProcessPoolExecutor() as executor:
            future_left = executor.submit(fdm_iteration_left, left_nx, ny, left_boundaries, v, left_rhs, xh, yh,
                                          fdm_iter)
            future_right = executor.submit(fdm_iteration_right, right_nx, ny, right_boundaries, v, right_rhs, xh, yh,
                                           fdm_iter)
            y_left = future_left.result()
            y_right = future_right.result()
        #y_left = fdm_iteration_left(left_nx, ny, left_boundaries, v, left_rhs, xh, yh, fdm_iter)
        #y_right = fdm_iteration_right(right_nx, ny, right_boundaries, v, right_rhs, xh, yh, fdm_iter)
        v[1] = y_left[left_x_mid_point]
        v[3] = y_right[right_x_mid_point]
        y = np.concatenate((y_left[:-1, :], y_right), axis=0)
        max_absolute_error = np.abs(exact_solution - y).max()
        errors.append(max_absolute_error)
        iteration_errors.append((i + 1, max_absolute_error))
        error = np.abs(exact_solution - y)
        print(f"iteration: {i + 1}, maximum absolute error:{max_absolute_error}")
        #update_plot_3d(left_x_nodes, right_x_nodes, y_nodes, exact_solution, y, i, xh, yh)
        plot_surface(left_x_nodes, right_x_nodes, y_nodes, exact_solution, y, i, xh, yh)
        #plot_error(left_x_nodes, right_x_nodes, y_nodes, error, i, xh, yh)
    #plot_iteration_errors(iteration_errors, xh, yh)

    min_error = min(errors)
    min_error_iteration = errors.index(min_error) + 1
    print(f"absolute error was minimal {min_error} at iteration {min_error_iteration}, ")


if __name__ == "__main__":
    pi = np.pi
    sin = np.sin
    cos = np.cos


    def u1(x, y):
        return x * y * sin(pi * y) * cos(pi * x)


    def f1(x, y):
        return 2 * x * pi * cos(pi * x) * cos(pi * y) - 2 * pi * y * sin(pi * y) * (sin(pi * x) + pi * x * cos(pi * x))


    def u2(x, y):
        return (1 - x) * y * sin(pi * y) * cos(pi * (1 - x))


    def f2(x, y):
        return 2 * pi * (x - 1) * cos(pi * x) * (cos(pi * y) - pi * y * sin(pi * y)) - 2 * pi * y * sin(pi * x) * sin(
            pi * y)


    initial_x_nodes = np.arange(0, 1.25 - 0.0001, 0.25)
    initial_y_nodes = np.arange(0, 1.5 - 0.0001, 1)

    h_x = 0.1
    h_y = 0.1

    finite_dif_iter = 5000
    iter = 10

    start = time.time()

    sss = non_local_contact_problem_iter(
        initial_x_nodes, initial_y_nodes,
        u1, u2, f1, f2, iter, finite_dif_iter, h_x, h_y
    )

    end = time.time()

    print(abs(end-start))

'''
    def u1(x, y):
        return 4 * x ** 4 + 8 * y ** 4

    def f1(x, y):
        return 48 * (x ** 2 + 2 * y ** 2)


    def u2(x1, y1):
        return 2 * x1 ** 3 + 8 * y1 ** 4


    def f2(x1, y1):
        return 12 * (x1 + 8 * y1 ** 2)'''
