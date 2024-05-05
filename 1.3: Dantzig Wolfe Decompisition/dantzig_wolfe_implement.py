import gurobipy as gp
from gurobipy import GRB
import random
from dantzig_wolfe import *
from gurobipy import Model, GRB


easy_constraints = easy_constraints_list('/Users/nedimozden/Desktop/CMOR442/HW4/constraints.txt')
file_path = '/Users/nedimozden/Desktop/CMOR442/HW4/22433.mps'
q_model = create_q(file_path, easy_constraints)

extreme_points = []
extreme_rays = []
for num in range(10):
    solve_q_with_rand_obj(q_model, extreme_points, extreme_rays)

coef = get_coefficients(file_path)

model, rhs_expr_pairs = create_rhs_expr_pairs(file_path, easy_constraints)

model, pi, sigma = create_dual(model, rhs_expr_pairs, extreme_points, extreme_rays, coef)

print(model.ObjVal)

q = create_q(file_path, easy_constraints)

pricing_answer, optimal_vars, is_bounded = pricing_problem(q, sigma, pi, coef, rhs_expr_pairs)

finished = evaluate_pricing_problem(pricing_answer, optimal_vars, is_bounded, extreme_points, extreme_rays)

while not finished:
    coef = get_coefficients(file_path)

    model, rhs_expr_pairs = create_rhs_expr_pairs(file_path, easy_constraints)

    model, pi, sigma = create_dual(model, rhs_expr_pairs, extreme_points, extreme_rays, coef)

    q = create_q(file_path, easy_constraints)

    pricing_answer, optimal_vars, is_bounded = pricing_problem(q, sigma, pi, coef, rhs_expr_pairs)

    finished = evaluate_pricing_problem(pricing_answer, optimal_vars, is_bounded, extreme_points, extreme_rays)

with open("HW4/optimization_result.txt", "w") as file:
    file.write(f"This optimal objective to the Dantzig-Wolfe Decomposition is: {model.ObjVal}\n")

