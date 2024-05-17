"""
best graphics:
import plotly.graph_objs as go
import plotly.io as pio
from plotly.subplots import make_subplots


def update_plot(nodes, exact_solution, approximations):
    fig = make_subplots(rows=1, cols=1, shared_xaxes=True, shared_yaxes=True)
    fig.add_trace(go.Scatter(x=nodes, y=exact_solution, mode='lines', name='Exact Solution'), row=1, col=1)
    for i, approx in enumerate(approximations):
        fig.add_trace(go.Scatter(x=nodes, y=approx, mode='lines', name=f'Approximation at iteration {i + 1}'),
                      row=1, col=1)
    fig.update_layout(title='Non-local Contact Problem Approximation', xaxis_title='x', yaxis_title='y')
    pio.show(fig)


this in iteration
approximations.append(y)

or after
update_plot(full_nodes, exact_solution, approximations)


"""


"""
import plotly.graph_objs as go
import plotly.io as pio

def plot_solution(nodes, exact_sol, approx):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=nodes, y=exact_sol, mode='lines', name='Exact Solution'))
    for i, k in enumerate(approx):
        fig.add_trace(go.Scatter(x=nodes, y=approx, mode='lines',
                                 name=f"Approximation at iteration {i + 1}"))
    fig.update_layout(title='Non-local Contact Problem Approximation', xaxis_title='x', yaxis_title='u(x)')
    pio.show(fig)

and put this in iteration:
        plot_solution(full_nodes, exact_solution, y)


"""

""""
    fig, ax = plt.subplots()
    ax.set_xlabel('x')
    ax.set_ylabel('u(x)')
    ax.set_title('Non-local Contact Problem Approximation')
    ax.grid(True)
    define this before iteration
    
    and
in non_local_problem define this after iteration
    def animate(frame):
        ax.clear()  # Clear the current plot
        ax.plot(full_nodes, exact_solution, label='Exact Solution')
        ax.plot(full_nodes, approximations[frame], label=f'Approximation at iteration {frame + 1}', alpha=0.7)
        ax.set_xlabel('x')
        ax.set_ylabel('u(x)')
        ax.set_title(f'Approximation at iteration {frame + 1}')
        ax.legend()
        ax.grid(True)

    ani = FuncAnimation(fig, animate, frames=iteration, interval=400)
    plt.gcf().set_dpi(100)
    plt.show()
"""



""" with open('iteration_errors.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Iteration', 'Max Absolute Error'])
        #for i, error in iteration_errors:
        #    writer.writerow([i, error])
        writer.writerow(iteration_errors[0])
        writer.writerow(iteration_errors[1000])
        writer.writerow((['Sufficient iteration', 'Min Absolute Error']))
        writer.writerow([min_error_iteration, min_error])"""