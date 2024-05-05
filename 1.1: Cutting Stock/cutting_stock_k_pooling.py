import gurobipy as gp
from gurobipy import GRB
import random
import os
import matplotlib.pyplot as plt

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


def k_pricing_problem(dual_variables, patterns, num_of_finals, raw_size, demands, k):
    knapsack_problem = gp.Model("knapsack_solution_pools")

    pi = {i: knapsack_problem.addVar(name=f"pi_{i}", lb=0, vtype=GRB.INTEGER) for i in range(num_of_finals)}

    # Objective function: maximize the sum of dual variables times the right-hand side of the primal constraints
    knapsack_problem.setObjective(gp.quicksum(dual_variables[i] * pi[i] for i in range(num_of_finals)), GRB.MAXIMIZE)

    # Constraints: one for each pattern in the primal problem
    for idx, pattern in enumerate(patterns):
        constraint_name = f"Ptrn_Cnst_{idx}"
        knapsack_problem.addConstr(gp.quicksum(demands[i][0] * pi[i] for i in range(num_of_finals)) <= raw_size, constraint_name)


    # Solution pool setting
    knapsack_problem.setParam('PoolSearchMode', 2)

    # Adding multiple columns at one iteration
    knapsack_problem.setParam('PoolSolutions', k)

    # Solve the model
    knapsack_problem.optimize()

    solutions = []

    for i in range(k):
        knapsack_problem.setParam('SolutionNumber', i)
        solutions.append([knapsack_problem.objVal, knapsack_problem.getAttr("Xn", knapsack_problem.getVars())])

    #print(solutions)



    return knapsack_problem, [solutions, knapsack_problem.Runtime, knapsack_problem.getAttr(gp.GRB.Attr.IterCount)]


def plot_comparisons(list_of_values, list_of_files, name_of_plot, attribute):
    x_values1, y_values1 = zip(*sorted(list_of_values[0].items()))
    x_values2, y_values2 = zip(*sorted(list_of_values[1].items()))
    x_values3, y_values3 = zip(*sorted(list_of_values[2].items()))

    plt.plot(x_values1, y_values1, label= list_of_files[0])

    # Plotting the second set of values
    plt.plot(x_values2, y_values2, label= list_of_files[1])

    # Plotting the third set of values
    plt.plot(x_values3, y_values3, label= list_of_files[2])

    # Adding labels and title
    plt.xlabel('K-Value')
    plt.ylabel(attribute)
    plt.title(name_of_plot)

    plt.savefig(name_of_plot + '.png')

    # Adding a legend
    plt.legend()


def main():
    with open('hw_1_pooling_sol.txt', 'w') as output_file:
        folder_path = "/Users/nedimozden/Desktop/CMOR442/HW1/"
        files = ["Scholl_CSP/N4W2B2R0.txt", "Scholl_CSP/N4W2B2R1.txt", "Scholl_CSP/N4W2B2R2.txt"]
        list_runtime = []
        list_iterations = []
        for file_name in files:
            # Iterate over each file and read its contents
            file_path = os.path.join(folder_path, file_name)
            dict_runtime = {}
            dict_iterations = {}
            list_runtime.append(dict_runtime)
            list_iterations.append(dict_iterations)
            for k in range(1, 20):
                total_runtime = 0
                total_iter = 0
                num_of_finals, raw_size, demands = read_instance(file_path)
                patterns = base_set(demands, raw_size, num_of_finals)
                
                #output_file.write("The base set  is: " + str(patterns) + "\n")
                primal_model = primalRMP(num_of_finals, patterns, demands)
                primal_model.update()
                primal_model.optimize()

                dual_model = dualRMP(patterns, demands, raw_size, num_of_finals)
                dual_model.update()
                # Solve the dual problem
                dual_model.optimize()

                dual_variables = [dual_model.getVarByName(f"y_{i}").X for i in range(num_of_finals)]
                #output_file.write("The dual variables are: " + str(dual_variables) + "\n")
                knapsack_model, knapsack_list = k_pricing_problem(dual_variables, patterns, num_of_finals, raw_size, demands, k)
                solutions, runtime, itercount = knapsack_list[0], knapsack_list[1], knapsack_list[2]
                total_runtime += runtime
                total_iter += itercount

                #output_file.write("The first pattern added is: " + str(pi_values) + "\n")
                knapsack_objective = knapsack_model.ObjVal
                #continuous column generation
                while (knapsack_objective > 1.0001):
                    for solution in solutions:
                        #output_file.write(str(solution) + "\n")
                        patterns.append(solution[1])
                    #run dualRMP
                    dual_model = dualRMP(patterns, demands, raw_size, num_of_finals)
                    dual_model.update()
                    dual_model.optimize()
                    dual_variables = [dual_model.getVarByName(f"y_{i}").X for i in range(num_of_finals)]

                    #run knapsack model
                    knapsack_model, knapsack_list = k_pricing_problem(dual_variables, patterns, num_of_finals, raw_size, demands, k)
                    solutions, runtime, itercount = knapsack_list[0], knapsack_list[1], knapsack_list[2]
                    total_runtime += runtime
                    total_iter += itercount
                    knapsack_objective = knapsack_model.ObjVal
                    

                primal_model = primalRMP(num_of_finals, patterns, demands)
                primal_model.update()
                primal_model.optimize()

                output_file.write("File " + "'" + str(file_name) + "' with k:" + str(k) + " columns added" + "\n")
                output_file.write("The total runtime is: " + str(total_runtime) + "\n")
                output_file.write("The total number of iterations are: " + str(total_iter) + "\n")
                dict_runtime[k] = total_runtime
                dict_iterations[k] = total_iter
                output_file.write("\n")
                output_file.flush()

        output_file.write(str(list_iterations) + "\n")
        output_file.write(str(list_runtime) + "\n")
        output_file.write("All instances have been solved")        
        output_file.close()

    # Extract x and y values from dictionaries
    plot_comparisons(list_runtime, files, "Runtime_Comparisons_of_Different_K_Values", "Runtime")

    plot_comparisons(list_iterations, files, "Iterations_Comparisons_of_Different_K_Values", "Iterations")

main()