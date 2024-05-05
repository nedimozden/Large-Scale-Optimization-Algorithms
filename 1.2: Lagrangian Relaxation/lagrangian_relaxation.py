import gurobipy as gp
from gurobipy import *

#read mps file
def read_mps_file(file_path):
    model = gp.read(file_path)
    counter = 1
    for constr in model.getConstrs():
        constr.ConstrName = str(counter)
        counter += 1
    model.update()
    return model
#read in hard constraints text file
def hard_constraints_list(file_path):
    hard_constraints = []
    with open(file_path, 'r') as file:
        for line in file:
            hard_constraints.append((line.strip()))
    return hard_constraints
#combination of changing to maximization and changing the constraints
def convert_to_proper_form(model, hc_list):
    changed_to_max = convert_to_maximization(model)
    #hc_list = change_constraints(model, hc_list)
    return model, hc_list, changed_to_max
#turn model into maximization problem if not already
def convert_to_maximization(model):
    if model.modelSense == GRB.MINIMIZE:
        changed_to_max = -1
        # Get the current objective function
        obj = model.getObjective()
        
        # Multiply the objective function by -1 to switch its direction
        new_obj = -1 * obj
        
        # Set the new objective as the model's objective
        model.setObjective(new_obj)
        model.modelSense = GRB.MAXIMIZE
    else:
        changed_to_max = 1
    
    model.update()
    return changed_to_max
#change all > constraints to <
def change_constraints(newmodel, list_hc):
    greater_than_constraints_to_modify = []

    # Step 1: Identify greater than constraints
    for constr in newmodel.getConstrs():
        constr_name = constr.ConstrName
        if constr.Sense == gp.GRB.GREATER_EQUAL:
            lhs_expr = newmodel.getRow(constr)
            rhs = constr.RHS
            
            # Prepare the negated version of the constraint
            negated_lhs_expr = -1 * lhs_expr
            negated_rhs = -1 * rhs

            # Store this constraint for modification
            greater_than_constraints_to_modify.append((constr, negated_lhs_expr, negated_rhs, constr_name))

    # Step 2: Modify the constraints
    for constr, negated_lhs_expr, negated_rhs, constr_name in greater_than_constraints_to_modify:
        # Remove the original constraint
        newmodel.remove(constr)

        # Add a new constraint with the negated LHS and RHS to make it a less than constraint
        newmodel.addConstr(negated_lhs_expr < negated_rhs, name = constr_name)
    # Don't forget to update the model after making changes
    newmodel.update()
    return list_hc        
#create sets of hard and easy constraints
def create_hard_and_easy_constraints(newmodel, hc_list):
    hard_constraints = []
    rhs_expr_pairs = []
    # Iterate through all constraints in the model
    for constr in newmodel.getConstrs():
        # Get the linear expression (LinExpr) associated with the constraint
        expr = newmodel.getRow(constr)
        constr_name = constr.ConstrName
        is_less_than = constr.Sense
        if constr_name in hc_list:
            rhs = constr.RHS
            rhs_expr_pairs.append([rhs, expr, is_less_than])
            hard_constraints.append(constr)
        
    for hardconstr in hard_constraints:
        #print("this is a hard constraint:", hardconstr.ConstrName)
        newmodel.remove(hardconstr)     
        newmodel.update()
    return newmodel, rhs_expr_pairs
#update the objective with the new lamda value
def calculate_zlr_optimal_x(newmodel, rhs_expr_pairs, lamda_val):
    objective = gp.LinExpr()
    objective.add(newmodel.getObjective())
    for i in range(len(rhs_expr_pairs)):
        objective.add(lamda_val[i] * (rhs_expr_pairs[i][0] - rhs_expr_pairs[i][1]))
    newmodel.setObjective(objective, sense=gp.GRB.MAXIMIZE)
    newmodel.update()
    
    newmodel.optimize()
    
    # Check if the model has an optimal solution
    if newmodel.status == gp.GRB.OPTIMAL:
        # Retrieve and return the values of the decision variables at the optimal solution
        optimal_x_values = [v.X for v in newmodel.getVars()]
        return newmodel, optimal_x_values, newmodel.ObjVal
    else:
        optimal_x_values = [v.X for v in newmodel.getVars()]
        # Handle cases where no optimal solution is found
        print("No optimal solution found.")
        return None
#compute lhs with optimal x's
def evaluate_expr(expr):
    # Initialize the value of the expression
    value = 0.0
    # Iterate over the terms in the expression
    for i in range(expr.size()):
        var = expr.getVar(i)  # Get the variable at position i in the expression
        coeff = expr.getCoeff(i)  # Get the coefficient at position i in the expression
        # Directly access the variable's value using the .X attribute
        value += coeff * var.X
    return value
#computes the subgradient at lambda!
def compute_subgradient(rhs_expr_pairs):
    subgradient = [0] * len(rhs_expr_pairs)
    zerovector = True
    for i in range(len(rhs_expr_pairs)):
        lhs = evaluate_expr(rhs_expr_pairs[i][1])
        subgradient[i] = rhs_expr_pairs[i][0] - lhs
        if subgradient[i] != 0:
            zerovector = False
    return subgradient, zerovector
#computes new lambda
def new_lambda(lambda_val, step_size, subgradient, rhs_expr_pairs):
    updated_lambda = lambda_val.copy()
    for i in range(len(lambda_val)):
        updated_lambda[i] = lambda_val[i] - (step_size * subgradient[i])
        if rhs_expr_pairs[i][2] == '<' and updated_lambda[i] < 0:
            updated_lambda[i] = 0
        if rhs_expr_pairs[i][2] == '>' and updated_lambda[i] > 0:
            updated_lambda[i] = 0
    return updated_lambda
#creates the zlr base objective
def create_base_zlr(file_path):
    #my new idea
    #read in old model, and make proper form
    og_model = read_mps_file(file_path)
    hc_list = hard_constraints_list("/Users/nedimozden/Desktop/CMOR442/HW3/constraints.txt")
    variable_model, hc_list, changed_to_max = convert_to_proper_form(og_model, hc_list) 
    variable_model, rhs_expr_pairs = create_hard_and_easy_constraints(variable_model, hc_list)
    variable_model.update()
    return variable_model, rhs_expr_pairs, changed_to_max
#tradition 1/t stepsize
def divergent_step_size(iteration):
    stepsize = 1.0 / iteration
    return stepsize
#calculates updated p_t value
def p_value(p_t, solution_list, iteration):
    #if the current obj = previous obj for the past 7 iterations -> p_t = p_(t-1) / 2
    if iteration >= 6 and (solution_list[iteration] == solution_list[iteration - 6]):
        p_t /= 2
    return p_t
#calculates the euclidian norm of the subgradient (used for nesterov step size)
def euclidian_norm_subgradient(subgradient):
    value = 0.0
    for component in subgradient:
        value += (component) ** 2
    return value
#calculates new nesterov step size
def nesterov_step_size(p_t, optimal, current, subgradient, changed_to_max):
    numerator = (current - (changed_to_max * optimal)) * p_t
    denominator = euclidian_norm_subgradient(subgradient)
    step_size = numerator / denominator
    return step_size

def calculate_sol(file_path):
    model = gp.read(file_path)
    model.update()
    model.optimize()
    return model.ObjVal
#puts everything together
def subgradient_algorithm(iterations, file_path, step_size_function, initial_lagrange):
    solution_list = [0] * iterations
    newmodel, rhs_expr_pairs, changed_to_max = create_base_zlr(file_path)
    ## Calculate ZLR with lambda vector of 1's
    ## Retrieve optimal X solution
    #lmbda = [initial_lagrange] * len(rhs_expr_pairs)
    lmbda = initial_lagrange
    new_model_copy, zlr_optimal, objVal = calculate_zlr_optimal_x(newmodel, rhs_expr_pairs, lmbda)
    best_sol = objVal
    solution_list[0] = objVal
    ## Use X values to calculate subgradient
    ## check if 0 vector is a subgradient
    subgradient, zero_vector = compute_subgradient(rhs_expr_pairs)
    # What I Need to Do
    ## Create a for loop that breaks when subgradient = 0 (maybe make True or False)
    optimal_feasible = calculate_sol(file_path)
    p_t = 2
    stepsize_dict = [0] * iterations
    for t in range(1, iterations, 1):
        if zero_vector:
            return (changed_to_max * best_sol)
        if step_size_function == divergent_step_size:
            stepsize = divergent_step_size(t)
        elif step_size_function == nesterov_step_size:
            p_t = p_value(p_t, solution_list, t)
            stepsize = nesterov_step_size(p_t, optimal_feasible, best_sol, subgradient, changed_to_max)
        stepsize_dict[t - 1] = stepsize
        lmbda = new_lambda(lmbda, stepsize, subgradient, rhs_expr_pairs)
        ## recalculate ZLR with new lambda
        newmodel, rhs_expr_pairs, changed_to_max = create_base_zlr(file_path)
        new_model_copy, zlr_optimal, compare_objVal = calculate_zlr_optimal_x(newmodel, rhs_expr_pairs, lmbda)
        solution_list[t] = compare_objVal
        ## check which solution is smaller
        ## if old solution is better, calculate new lambda by changing step size
        ## if new solution is better, calculate new subgradient
        subgradient, zero_vector = compute_subgradient(rhs_expr_pairs)
        if t == 1:
            fs = subgradient
        if compare_objVal < best_sol:
            best_sol = compare_objVal
    return (changed_to_max * best_sol), solution_list, stepsize_dict, fs

final_answer, solutionlist, stepsizedict, fs = subgradient_algorithm(20, "/Users/nedimozden/Desktop/CMOR442/HW3/mps_files/22433.mps", nesterov_step_size, [1, -1, 1, -1, 1])
real_answer = calculate_sol("/Users/nedimozden/Desktop/CMOR442/HW3/mps_files/22433.mps")
print(str(final_answer))
print(str(real_answer))
print(len(stepsizedict))