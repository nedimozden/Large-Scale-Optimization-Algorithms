import gurobipy as gp
from gurobipy import GRB
import random

def easy_constraints_list(file_path):
    easy_constraints = []
    with open(file_path, 'r') as file:
        for line in file:
            easy_constraints.append((line.strip()))
    return easy_constraints

def rename_constraints(file_path):
    model = gp.read(file_path)
    counter = 1
    for constr in model.getConstrs():
        constr.ConstrName = str(counter)
        counter += 1
    model.update()
    varCounter = 0
    for var in model.getVars():
        var.VarName = str(varCounter)
        varCounter +=1
    model.update()
    return model

#read mps file
def create_q(file_path, easy_constraints): 
    model = rename_constraints(file_path)
    all_constraints = model.getConstrs()
    for constr in all_constraints:
        if constr.ConstrName in easy_constraints:
            model.remove(constr)

    model.update()
    return model

def create_mp(file_path, easy_constraints): 
    model = rename_constraints(file_path)
    all_constraints = model.getConstrs()
    for constr in all_constraints:
        if constr.ConstrName in easy_constraints:
            model.remove(constr)

    model.update()
    return model

def solve_q_with_rand_obj(model, extreme_points, extreme_rays):
    """
    Load a model, modify it by keeping only specified hard constraints,
    generate a random objective function, and solve the model.

    Parameters:
    - model: a Gurobi model instance.
    - extreme_points: list, stores the lists of optimal x values if the objective is bounded.
    - extreme_rays: list, stores the lists of x values if the objective is unbounded.

    Returns:
    - A list of results, where each entry indicates the status of the model.
    """

    # List to store the status of each iteration
    results = []

    variables = model.getVars()

    # Generate a random objective
    coefficients = [random.uniform(-10, 10) for _ in variables]
    objective = gp.LinExpr(coefficients, variables)
    model.setObjective(objective, GRB.MAXIMIZE)

    # Optimize the model
    model.optimize()

    # Store the result as a list of variable values
    if model.status == GRB.OPTIMAL:
        # Collecting optimal values of variables
        result = [v.X for v in variables]
        if result in extreme_points:
            results.append('Optimal solution added to extreme points.')
            return results
        else:
            extreme_points.append(result)
            results.append('Optimal solution added to extreme points.')
    elif model.status == GRB.UNBOUNDED:
        # Collecting values contributing to the unboundedness (typically coefficients of the unbounded ray)
        result = [v.UnbdRay for v in variables]
        if result in extreme_rays:
            results.append('Unbounded solution added to extreme rays.')
            return results
        else:
            extreme_rays.append(result)
            results.append('Unbounded solution added to extreme rays.')
    else:
        result = []
        results.append('Model not optimal or infeasible.')

    # Reset the model to change objective in the next iteration without re-loading
    model.reset()
    return model

def get_coefficients(file_path):
    model = gp.read(file_path)
    c = [var.getAttr("obj") for var in model.getVars()]
    return c

def create_rhs_expr_pairs(file_path, easy_constraints):
    model = rename_constraints(file_path)

    hard_constraints = []
    rhs_expr_pairs = []
    for constr in model.getConstrs():
        if constr.ConstrName in easy_constraints:

            model.remove(constr)
            model.update()
            # Get the linear expression (LinExpr) associated with the constraint
            
        else:
            hard_constraints.append(constr)
            expr = model.getRow(constr)
            rhs = constr.RHS  # Right-hand side of the constraint
            a_1 = [0] * model.NumVars
            for i in range(expr.size()):
                var = expr.getVar(i)  # Get the variable at position i in the expression
                name = var.VarName
                coeff = expr.getCoeff(i)
                a_1[int(name)] = coeff
            rhs_expr_pairs.append((rhs, a_1))
    for hard_constr in hard_constraints:
        model.remove(hard_constr)
        model.update()

    return model, rhs_expr_pairs

def create_dual(model, rhs_expr_pairs, extreme_points, extreme_rays, coef):
    # Create dual variables
    pi = [model.addVar(lb=0, vtype=GRB.CONTINUOUS, name=f"pi_{i}") for i in range(len(coef))]
    sigma = model.addVar(vtype=GRB.CONTINUOUS, name="sigma")
    
    # Update model to integrate new variables
    model.update()

    # Set the objective function
    objective = gp.quicksum(pi[i] * rhs_expr_pairs[i][0] for i in range(len(rhs_expr_pairs))) + sigma
    model.setObjective(objective, GRB.MINIMIZE)

    # Add constraints for each extreme point
    for point in extreme_points + extreme_rays:
        for rhs, a_1 in rhs_expr_pairs:
            constraint_expr_lhs = gp.quicksum(pi[j] * a_1[j] * point[j] for j in range(len(a_1)))
            constraint_expr_rhs = gp.quicksum(point[j] * coef[j] for j in range(len(coef)))
            model.addConstr(constraint_expr_lhs + sigma >= constraint_expr_rhs)

    # Optimize the model
    model.optimize()

    # Check the solution status and return results accordingly
    if model.status == GRB.OPTIMAL:
        pi_values = [var.X for var in pi]
        sigma_value = sigma.X
        return model, pi_values, sigma_value
    else:
        print(f"No optimal solution found. Model status: {model.status}")
        return None



from gurobipy import GRB, Model

def create_dual_model(primal_model):
    # Create a new model for the dual
    dual_model = Model("dual")

    # Create dual variables for each constraint in the primal
    dual_vars = {}
    for constr in primal_model.getConstrs():
        if constr.Sense == GRB.LESS_EQUAL:
            # For 'less than or equal to' constraints in primal, dual variable should be non-negative
            dual_vars[constr.ConstrName] = dual_model.addVar(lb=0, obj=constr.RHS, name=f"dual_{constr.ConstrName}")
        elif constr.Sense == GRB.GREATER_EQUAL:
            # For 'greater than or equal to' constraints in primal, dual variable should be non-positive
            dual_vars[constr.ConstrName] = dual_model.addVar(ub=0, obj=constr.RHS, name=f"dual_{constr.ConstrName}")
        else:
            # For equality constraints in primal, dual variable is unrestricted
            dual_vars[constr.ConstrName] = dual_model.addVar(obj=constr.RHS, name=f"dual_{constr.ConstrName}")

    # Add constraints to the dual model based on the primal variables
    primal_vars = primal_model.getVars()
    for var in primal_vars:
        # Collect the dual variables and coefficients associated with the primal variable 'var'
        dual_expr = gp.quicksum(dual_vars[constr.ConstrName] * constr.getConstr(primal_model, constr.index, var.VarName)
                             for constr in primal_model.getConstrs()
                             if var.VarName in constr.getAttr('VarName'))

        # Add the dual constraint, direction depends on whether primal is minimization or maximization
        if primal_model.ModelSense == GRB.MINIMIZE:
            dual_model.addConstr(dual_expr <= var.Obj, name=f"dual_constr_{var.VarName}")
        else:
            dual_model.addConstr(dual_expr >= var.Obj, name=f"dual_constr_{var.VarName}")

    dual_model.update()
    return dual_model


def pricing_problem(q, sigma_optimal, pi_optimal, coef, rhs_expr_pairs):
    new_objective = gp.LinExpr()
    vars = q.getVars()
    for rhs, a_1 in (rhs_expr_pairs):
        new_objective += gp.quicksum((coef[j] - (pi_optimal[j] * a_1[j]) * vars[j]) for j in range(len(a_1)))
    q.setObjective(new_objective, gp.GRB.MAXIMIZE)

    q.update()
    q.optimize()

    
    # Initialize the is_bounded flag
    is_bounded = True

    # Check the optimization status
    if q.status == gp.GRB.OPTIMAL:
        pricing_answer = q.objVal - sigma_optimal
        # Retrieve and store the optimal values of the variables
        optimal_vars = [var.x for var in vars]
    elif q.status == gp.GRB.UNBOUNDED:
        # Set is_bounded to False if the model is unbounded
        is_bounded = False
        pricing_answer = float('inf')  # Or some other indicative value for unbounded
        optimal_vars = [v.UnbdRay for v in vars]
    else:
        raise Exception("Optimization was unsuccessful. Status code:", q.status)
    
    return pricing_answer, optimal_vars, is_bounded


def evaluate_pricing_problem(pricing_answer, optimal_vars, is_bounded, extreme_points, extreme_rays):
    finished = False
    if not is_bounded:
        extreme_rays.append(optimal_vars)
        return finished
    elif pricing_answer > 0 and is_bounded:
        extreme_points.append(optimal_vars)
        return finished
    else:
        finished = True
        return finished