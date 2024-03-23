""""import numpy as np

nodes = np.array([0, 1/4, 1/2, 3/4, 1])
a = 0.25
b = 0.5
h = nodes[1] - nodes[0]
def exact_u1(x):
    return 1/48*(8 * x**4 + 71 * x)


def exact_u2(x):
    return 1/6 * (2 * x**3 + 15 * x**2 - 35 * x + 18)


def f_1(x):
    return 2 * x**2


def f_2(x):
    return 2 * x + 5


f_0 = exact_u1(nodes[2]) - (a * exact_u1(nodes[1]) + b * exact_u2(nodes[3]))


errors = {}

def iterate(u1, u2, k):
    for i in range(k):
        u1_new_half = a * u1 + b * u2 + f_0
        u2_new_half = a * u1 + b * u2 + f_0
        u1_new_quarter = (exact_u1(nodes[0]) + u1_new_half - f_1(nodes[1]) * h**2) * 0.5
        u2_new_three_quarters = (exact_u2(nodes[4]) + u2_new_half - f_2(nodes[3]) * h**2) * 0.5

        u1 = u1_new_quarter
        u2 = u2_new_three_quarters

        u_half = a * u1 + b * u2 + f_0 # Using the value at x = 1/2
        exact_result = exact_u1(nodes[2])
        errors[i+1] = abs(exact_result - u_half)
        print(f"Iteration {i + 1}: u(1/2) = {u_half}")
        print(f"Absolute Error {i + 1}: {errors[i+1]}")
    iteration_minimum_error = min(errors, key=errors.get)
    print(f"iteration : {iteration_minimum_error} : we get minimum error: {min(errors.values())}")
    print(f"exact value in 1/2: {exact_u1(nodes[2])}")


# Initial values and parameters
u1_initial = 0
u2_initial = 0
iteration = 7  # number of iterations

# Perform iterations
iterate(u1_initial, u2_initial, iteration)
""""