import gurobipy as gp
from gurobipy import GRB
import os

def read_instance(file_path):
    num_of_finals = 0
    raw_size = 0
    demands = []

    with open(file_path, 'r') as file:
        # this is getting the first number on the text file and setting it equal to the number of finals
        num_of_finals = int(file.readline().strip())  
        # this is getting the second number on the text file and setting it euqal to the size of the rolls
        raw_size = int(file.readline().strip())  

        # for the remaining lines in the text file, we split it into the final_size needed and its demand
        for line in file:
            size_needed, demand_size = map(int, line.split())
            demands.append([size_needed, demand_size])

    return num_of_finals, raw_size, demands


def base_set(demands, raw_size, num_of_finals):
    #create the initial set of patterns
    patterns = []
    for index, pair in enumerate(demands):
        #max number of each final that can fit in one raw
        max_value = raw_size // pair[0]
        new_list = [0] * num_of_finals
        new_list[index] = max_value
        patterns.append(new_list)

    return patterns


def primalRMP(num_of_finals, patterns, demands):
    primalRMP = gp.Model("PrimalRMP")

    # Decision variables (one for each pattern in the primal problem)
    x = {j: primalRMP.addVar(name=f"x_{j}", lb=0, vtype=GRB.CONTINUOUS) for j in range(len(patterns))}

    # Objective function: minimize the sum of x variables
    primalRMP.setObjective(gp.quicksum(x[j] for j in range(num_of_finals)), sense=GRB.MINIMIZE)

    # Constraints: one for each demand in the primal problem
    for i in range(num_of_finals):
        constraint_expr = gp.quicksum(patterns[j][i] * x[j] for j in range(len(patterns))) >= demands[i][1]
        constraint_name = f"constraint_{i}"
        primalRMP.addConstr(constraint_expr, constraint_name)
        primalRMP.update()
        #print(f"Added constraint: {constraint_name}, Expression: {constraint_expr}")
        

    return primalRMP


def dualRMP(patterns, demands, raw_size, num_of_finals):
    # Create Gurobi model for the dual problem
    dualRMP = gp.Model("DualRMP")

    # Dual variables (one for each constraint in the primal problem)
    y = {i: dualRMP.addVar(name=f"y_{i}", lb=0, vtype=GRB.CONTINUOUS) for i in range(num_of_finals)}

    # Objective function: maximize the sum of dual variables times the right-hand side of the primal constraints
    dualRMP.setObjective(gp.quicksum(demands[i][1] * y[i] for i in range(num_of_finals)), GRB.MAXIMIZE)

    # Constraints: one for each pattern in the primal problem
    for idx, pattern in enumerate(patterns):
        constraint_name = f"Ptrn_Cnst_{idx}"
        dualRMP.addConstr(gp.quicksum(pattern[i] * y[i] for i in range(num_of_finals)) <= 1, constraint_name)


    return dualRMP


def knapsack_problem(dual_variables, patterns, num_of_finals, raw_size, demands):
    knapsack_problem = gp.Model("KnapsackProblem")

    # Dual variables (one for each constraint in the primal problem)
    pi = {i: knapsack_problem.addVar(name=f"pi_{i}", lb=0, vtype=GRB.INTEGER) for i in range(num_of_finals)}

    # Objective function: maximize the sum of dual variables times the right-hand side of the primal constraints
    knapsack_problem.setObjective(gp.quicksum(dual_variables[i] * pi[i] for i in range(num_of_finals)), GRB.MAXIMIZE)

    # Constraints: one for each pattern in the primal problem
    for idx, pattern in enumerate(patterns):
        constraint_name = f"Ptrn_Cnst_{idx}"
        knapsack_problem.addConstr(gp.quicksum(demands[i][0] * pi[i] for i in range(num_of_finals)) <= raw_size, constraint_name)

    return knapsack_problem
 

def main():

    with open('hw_1_optimalsolution_file.txt', 'w') as output_file:
        folder_path = "/Users/nedimozden/Desktop/CMOR442/HW1/Scholl_CSP"

        # Get a list of all files in the folder
        file_list = [f for f in os.listdir(folder_path) if f.endswith('.txt')]

        # Iterate over each file and read its contents
        for file_name in file_list:
            file_path = os.path.join(folder_path, file_name)
            num_of_finals, raw_size, demands = read_instance(file_path)
            patterns = base_set(demands, raw_size, num_of_finals)
            primal_model = primalRMP(num_of_finals, patterns, demands)
            primal_model.update()
            primal_model.optimize()
            dual_model = dualRMP(patterns, demands, raw_size, num_of_finals)
            dual_model.update()
            # Solve the dual problem
            dual_model.optimize()

            dual_variables = [dual_model.getVarByName(f"y_{i}").X for i in range(num_of_finals)]

            knapsack_model = knapsack_problem(dual_variables, patterns, num_of_finals, raw_size, demands)
            knapsack_model.update()
            knapsack_model.optimize()
            # Extract pi variable values
            pi_values = [abs(knapsack_model.getVarByName(f"pi_{i}").X) for i in range(num_of_finals)]

            #continuous column generation
            while (pi_values not in patterns):

                primal_model = primalRMP(num_of_finals, patterns, demands)
                primal_model.update()
                primal_model.optimize()

                patterns.append(pi_values)
                #run dualRMP
                dual_model = dualRMP(patterns, demands, raw_size, num_of_finals)
                dual_model.update()
                dual_model.optimize()
                dual_variables = [dual_model.getVarByName(f"y_{i}").X for i in range(num_of_finals)]

                #run knapsack model
                knapsack_model = knapsack_problem(dual_variables, patterns, num_of_finals, raw_size, demands)
                knapsack_model.update()
                knapsack_model.optimize()
                pi_values = [abs(knapsack_model.getVarByName(f"pi_{i}").X) for i in range(num_of_finals)]

            primal_model = primalRMP(num_of_finals, patterns, demands)
            primal_model.update()
            primal_model.optimize()

            x_variables = [primal_model.getVarByName(f"x_{j}").X for j in range(len(patterns))]

            output_file.write("The optimal solution for file " + "'" + str(file_name) + "' is: %.2f" % primal_model.objVal + "\n")
            """for j in range(len(patterns)):
                if x_variables[j] == 0:
                    continue
                else:
                    output_file.write("\t" + str(x_variables[j]) + " copies of pattern " + str(patterns[j]) + "\n")
                    output_file.write("\n")
                    """
            output_file.write("\n")
            output_file.flush()


    output_file.write("All instances have been solved")        
    output_file.close()

main()





        
        

    
