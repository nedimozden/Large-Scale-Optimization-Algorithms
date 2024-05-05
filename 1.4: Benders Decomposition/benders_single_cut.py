import gurobipy as gp
from gurobipy import GRB
import math

cases = [[3, 3.6, 24], [2.4, 3, 20], [2, 2.4, 16], [2.2, 3.4, 17], [2.8, 2.6, 19], [2.1, 2.5, 22]]
probabilities = [0.1, 0.2, 0.15, 0.35, 0.15, 0.05]
penalties = [-200, -240, 0, -6000]



def compute_dual_parameters(variables, cases):
    # Initialize a new optimization model for calculating recourse
    recourse_model = gp.Model("RecourseOptimization")
    
    # Define recourse variables to be optimized
    recourse_vars = {
        'w1': recourse_model.addVar(lb=0, name="w1", vtype=GRB.CONTINUOUS),
        'w2': recourse_model.addVar(lb=0, name="w2", vtype=GRB.CONTINUOUS),
        'w3': recourse_model.addVar(lb=0, name="w3", vtype=GRB.CONTINUOUS),
        'w4': recourse_model.addVar(lb=0, name="w4", vtype=GRB.CONTINUOUS)
    }
    
    # Set the objective function to minimize negative impacts from recourse variables
    recourse_model.setObjective(-170 * recourse_vars['w1'] - 150 * recourse_vars['w2'] - 36 * recourse_vars['w3'] - 10 * recourse_vars['w4'], GRB.MINIMIZE)
    
    # Define constraints based on the scenario specifics and decision variables
    recourse_model.addConstr(cases[0] * variables[0] - recourse_vars['w1'] >= 200, "Constraint1")
    recourse_model.addConstr(cases[1] * variables[1] - recourse_vars['w2'] >= 240, "Constraint2")
    recourse_model.addConstr(recourse_vars['w3'] + recourse_vars['w4'] - cases[2] * variables[2] <= 0, "Constraint3")
    recourse_model.addConstr(recourse_vars['w3'] <= 6000, "Constraint4")
    
    # Solve the optimization problem
    recourse_model.optimize()
    
    # Retrieve dual variables from the optimization model
    duals = recourse_model.getAttr(GRB.Attr.Pi)

    return recourse_model.ObjVal, duals

def calculate_optimality(vars):
    total_recourse = []
    expected_cost = 0
    weighted_coeffs = [0, 0, 0]
    
    for index, case in enumerate(cases):
        # Compute recourse and dual variables for each case
        recourse_value, dual_variables = compute_dual_parameters(vars, case)
        total_recourse.append(recourse_value)
        
        # Calculate expected cost contribution for this scenario
        scenario_contribution = probabilities[index] * sum(dual_variables[j] * penalties[j] for j in range(len(penalties)))
        expected_cost += scenario_contribution
        
        # Update coefficients for constraints in the optimization model
        for j in range(len(weighted_coeffs)):
            weighted_coeffs[j] += probabilities[index] * dual_variables[j] * case[j]

    # Calculate weighted sum of recourses and adjust by coefficients
    weighted_sum_recourses = sum(total_recourse[i] * probabilities[i] for i in range(len(total_recourse)))
    cut_value = weighted_sum_recourses - sum(weighted_coeffs)
    
    return expected_cost, weighted_coeffs, cut_value

def evaluate_feasibility(decision_vars, scenario):
    # Create a model to evaluate feasibility of the decisions given the scenario constraints
    feasibility_test = gp.Model("FeasibilityEvaluation")

    # Declare recourse and slack variables in a structured way
    recourse_vars = {f'w{i+1}': feasibility_test.addVar(lb=0, name=f"w{i+1}", vtype=GRB.CONTINUOUS) for i in range(4)}
    slack_vars_plus = {f'nuPlus{i+1}': feasibility_test.addVar(lb=0, name=f"nuPlus{i+1}", vtype=GRB.CONTINUOUS) for i in range(4)}
    slack_vars_minus = {f'nuMinus{i+1}': feasibility_test.addVar(lb=0, name=f"nuMinus{i+1}", vtype=GRB.CONTINUOUS) for i in range(4)}
    
    # Define the objective function to minimize the total slack (penalty for deviation from constraints)
    feasibility_test.setObjective(sum(slack_vars_plus.values()) + sum(slack_vars_minus.values()), GRB.MINIMIZE)

    # Set up constraints to check if the scenario can be met with the current decision variables
    feasibility_test.addConstr(slack_vars_plus['nuPlus1'] - slack_vars_minus['nuMinus1'] - scenario[0] * decision_vars[0] + recourse_vars['w1'] <= -200, "Constraint1")
    feasibility_test.addConstr(slack_vars_plus['nuPlus2'] - slack_vars_minus['nuMinus2'] - scenario[1] * decision_vars[1] + recourse_vars['w2'] <= -240, "Constraint2")
    feasibility_test.addConstr(slack_vars_plus['nuPlus3'] - slack_vars_minus['nuMinus3'] - scenario[2] * decision_vars[2] + recourse_vars['w3'] + recourse_vars['w4'] <= 0, "Constraint3")
    feasibility_test.addConstr(slack_vars_plus['nuPlus4'] - slack_vars_minus['nuMinus4'] + recourse_vars['w4'] <= 6000, "Constraint4")

    # Solve the feasibility problem to find minimal slack needed
    feasibility_test.optimize()

    # Get dual values to understand the marginal prices of constraints
    dual_values = feasibility_test.getAttr(GRB.Attr.Pi)

    return feasibility_test.ObjVal, dual_values

def create_master_problem():
    master_problem = gp.Model("L Shaped Method")
    decision_var_1 = master_problem.addVar(name="decision_var_1", lb=0, vtype=GRB.CONTINUOUS)
    decision_var_2 = master_problem.addVar(name="decision_var_2", lb=0, vtype=GRB.CONTINUOUS)
    decision_var_3 = master_problem.addVar(name="decision_var_3", lb=0, vtype=GRB.CONTINUOUS)
    master_problem.setObjective(150 * decision_var_1 + 230 * decision_var_2 + 260 * decision_var_3, GRB.MINIMIZE)
    master_problem.addConstr(decision_var_1 + decision_var_2 + decision_var_3 <= 500)
    master_problem.update()
    return master_problem, decision_var_1, decision_var_2, decision_var_3

def run_optimization(master_problem, var_1, var_2, var_3):
    iterations = 0
    while True:
        master_problem.optimize()
        current_values = [var_1.X, var_2.X, var_3.X]
        
        # Update constraints based on feasibility
        while True:
            total_penalty = 0
            for t in cases:
                feasibility, duals = evaluate_feasibility(current_values, t)
                coeffs = [t[j] * duals[j] for j in range(3)]
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
            theta_val = -math.inf
            theta = master_problem.addVar(name='theta', lb=-math.inf, vtype=GRB.CONTINUOUS)
            master_problem.setObjective(150 * var_1 + 230 * var_2 + 260 * var_3 + theta, GRB.MINIMIZE)
            master_problem.update()
        else:
            theta_val = theta.X

        # Calculate optimality and update constraints accordingly
        e, cut, w = calculate_optimality(current_values)

        if w < theta_val:
            for constr in master_problem.getConstrs():
                print(str(master_problem.getRow(constr)) + str(constr.Sense) + str(constr.RHS))
            return current_values, iterations
        master_problem.addConstr(cut[0] * var_1 + cut[1] * var_2 - cut[2] * var_3 + theta >= -e)
        print(e)
        master_problem.update()
        iterations += 1


# Run L-shaped method
master_problem, x1, x2, x3 = create_master_problem()
print(run_optimization(master_problem, x1, x2, x3))
