import gurobipy as gp
from gurobipy import GRB
import math
from benders_single_cut import *

cases = [[3, 3.6, 24], [2.4, 3, 20], [2, 2.4, 16], [2.2, 3.4, 17], [2.8, 2.6, 19], [2.1, 2.5, 22]]
probabilities = [0.1, 0.2, 0.15, 0.35, 0.15, 0.05]
penalties = [-200, -240, 0, -6000]


def calculate_optimality(x_values, index):
    case = cases[index]
    recourse_value, dual_variables = compute_dual_parameters(x_values, case)   
    expected_cost = probabilities[index] * sum(dual_variables[j] * penalties[j] for j in range(len(penalties)))

    weighted_coeffs = [probabilities[index] * dual_variables[0] * case[0],
                       probabilities[index] * dual_variables[1] * case[1],
                       probabilities[index] * dual_variables[2] * case[2]]

    w_value = probabilities[index] * recourse_value - weighted_coeffs[0] - weighted_coeffs[1] - weighted_coeffs[2]
    return expected_cost, weighted_coeffs, w_value

def run_optimization(master_problem, var_1, var_2, var_3):
    iterations = 0
    while iterations < 10:
        master_problem.optimize()
        current_values = [var_1.X, var_2.X, var_3.X]
        
        # Update constraints based on feasibility
        while True:
            total_penalty = 0
            for case in cases:
                feasibility, duals = evaluate_feasibility(current_values, case)
                coeffs = [case[j] * duals[j] for j in range(3)]
                total_penalty += feasibility
                if feasibility != 0:
                    master_problem.addConstr(-coeffs[0] * var_1 - coeffs[1] * var_2 - coeffs[2] * var_3 >= feasibility)
                    master_problem.update()
            if total_penalty > 0:
                master_problem.optimize()
                current_values = [var_1.X, var_2.X, var_3.X]
            else:
                break

        # Initialize or update theta variable
        if iterations == 0:
            theta_vals = [-float('inf') for _ in range(6)]
            theta_vars = [master_problem.addVar(name=f"theta_{i}", lb=-float('inf'), vtype=GRB.CONTINUOUS) for i in range(6)]
            master_problem.setObjective(150 * var_1 + 230 * var_2 + 260 * var_3 + gp.quicksum(theta_vars), GRB.MINIMIZE)
            master_problem.update()
        else:
            theta_vals = [theta.X for theta in theta_vars]

        x = 0
        for i, variable in enumerate(theta_vars):
            # Calculate optimality and update constraints accordingly
            e, cut, w = calculate_optimality(current_values, i)
            if w > theta_vals[i]:
                master_problem.addConstr(cut[0] * var_1 + cut[1] * var_2 - cut[2] * var_3 + variable >= -e)
                x += 1

        # if all of the ws are less than their respective theta, we are done
        if x == 0:
            return current_values, iterations

        master_problem.update()
        iterations += 1


# Run L-shaped method
master_problem, x1, x2, x3 = create_master_problem()
print(run_optimization(master_problem, x1, x2, x3))