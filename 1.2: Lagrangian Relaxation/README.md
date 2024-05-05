# Lagrangian Relaxation with Divergent and Nesterov Step Sizes

This repository contains implementations of Lagrangian Relaxation with two types of step sizes: divergent and Nesterov. Lagrangian Relaxation is a powerful technique used in optimization problems to obtain approximate solutions by relaxing constraints and solving simpler subproblems iteratively. The choice of step size influences the convergence and efficiency of the algorithm.

## Introduction

Lagrangian Relaxation is a technique commonly used to solve optimization problems with complex constraints. By relaxing some constraints and solving simplified subproblems iteratively, it can often provide near-optimal solutions efficiently. The choice of step size in Lagrangian Relaxation affects the convergence rate and the quality of the approximate solutions obtained.

In this repository, we provide implementations of Lagrangian Relaxation with two types of step sizes:

### Divergent Step Size: 
In this approach, the step size increases iteratively, allowing for faster convergence initially but risking overshooting the optimal solution.

### Nesterov Step Size: 
Nesterov accelerated gradient (NAG) is a variant of gradient descent that uses momentum to accelerate convergence. Applying Nesterov step size in Lagrangian Relaxation can improve convergence speed and stability.
