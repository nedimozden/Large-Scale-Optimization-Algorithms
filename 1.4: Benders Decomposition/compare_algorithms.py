import gurobipy as gp
from gurobipy import GRB
import time
from benders_single_cut import create_master_problem as create_single_master_problem, run_optimization as run_single_optimization
from benders_multi_cut import create_master_problem as create_multi_master_problem, run_optimization as run_multi_optimization

import time

# Single Cut Algorithm
single_cut_start_time = time.time()
master_problem_single, x1_single, x2_single, x3_single = create_single_master_problem()
solution_single, iterations_single = run_single_optimization(master_problem_single, x1_single, x2_single, x3_single)
single_cut_end_time = time.time()
single_cut_solution_time = single_cut_end_time - single_cut_start_time

# Multi Cut Algorithm
multi_cut_start_time = time.time()
master_problem_multi, x1_multi, x2_multi, x3_multi = create_multi_master_problem()
solution_multi, iterations_multi = run_multi_optimization(master_problem_multi, x1_multi, x2_multi, x3_multi)
multi_cut_end_time = time.time()
multi_cut_solution_time = multi_cut_end_time - multi_cut_start_time

# Writing results to a text file
with open("HW5/algorithm_comparison_results.txt", "w") as file:
    file.write("Algorithm Comparison\n\n")
    
    file.write("Single Cut Algorithm:\n")
    file.write(f"Solution: {solution_single}\n")
    file.write(f"Number of Iterations: {iterations_single}\n")
    file.write(f"Solution Time: {single_cut_solution_time} seconds\n\n")
    
    file.write("Multi Cut Algorithm:\n")
    file.write(f"Solution: {solution_multi}\n")
    file.write(f"Number of Iterations: {iterations_multi}\n")
    file.write(f"Solution Time: {multi_cut_solution_time} seconds\n")

