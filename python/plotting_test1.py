import sys

import numpy as np
import gmpy2 as gm

np.set_printoptions(threshold=sys.maxsize)

# Methods to find values of alpha, beta [boundary values]
from recursions import alpha, beta, gamma

# Methods for generating values and normal derivatives of polynomials
from monomials import p_jk, p_jk_addr_list, generate_W, generate_T, \
    generate_norm_W, generate_norm_T

# Methods for plotting functions on SG
from plotting import plot_monomial, plot_easy_basis, plot_op, \
    plot_general

# Methods for evaluating Chebyshev polynomials
from chebyshev import eval_poly_h1, plot_h1, extremal_val,  \
    extremal_val_h1, max_h1_val_family, max_h2_val_family, plot_hn, \
    eval_poly_hn, extremal_val_hn, max_h1_full_family


####### Outline of the testing functionalities here (use search!)
# Section I. Basic Plotting
# 1. Plotting the easy basis/monomial/SOP
# 2. Plotting a general polynomial

# Section II. Maximal, Minimal, Sup-norms
# 1. Code for finding sup-norm of general polynomials

# Section III. Finding Chebyshev Polynomials

# Section IV. Finding Optimality Conditions


###-------------------- Section I. Plotting -------------------###

########## 1. Code for Plotting the Easy Basis/Monomial/SOP ##########

### Parameters
# num represent the highest power to be plotted
# k represent the type of monomial {0,1,2} -> {1, 2, 3}
# num = 3
# k = 0

# # Example A. To plot the monomials, use plot_monomial
# plot_monomial(num, k)


# # Example B. To plot the easy basis, use plot_easy_basis
# plot_easy_basis(num, k)


# # Example C. This plots the Sobolov Orthogonal Polynomials
# plot_op(num, k)


########## 2. Code for Plotting a general polynomial in H1 ##########

# # Example D: Compute P_{13} + 1/2*P_{03}
# plot_h1(0, 0, 1/2, 0, 0, 1, 7)



###---------- Section II. Maximal, Minimal, Sup-norms ----------###

######## 1. Code for finding sup-norm of general polynomials ########

### Example A: Sup-norm of f_{j0}
### We can see from the following code that f_{j0} always take the maximum
### norm at point F_0F_1q_2=F_0F_2q_1

# # Initiaize the parameters
# j = 3
# num_points = 10

# # Generate the values of the polynomial
# W = generate_W(7, j+1, False)
# fj1_val = W[0, :, j]

# # Find the extremal value
# extremal_val(fj1_val, num_points, flag=3, level=7)


### Example B: Sup-norm of f_{j+1, k} / f_{j, k} 
### We can see that the sup-norm is achieved at F_{1}q_{2}

# # Initiaize the parameters
# j = 7
# num_points = 10

# # Generate the values of the polynomial
# W = generate_W(7, j+2, False)
# fj0_seq = W[0, :, j]
# fj1_0_seq = W[0, :, j+1]
# prop = np.abs(np.divide(fj1_0_seq, fj0_seq))

# # Find the extremal value
# extremal_val(prop, num_points, flag=3, level=7)


### Note that the maximum/minimum/smallest norm can be similarly 
### obtained if we change the flag variable in extreme_val.



###-------------- III. Finding Chebyshev Polynomials ---------------###

### Example 0: Chebyshev Polynomial of order 0
### For completeness, we also graph the Chebyshev polynomials on level 0
# # T_{01} = P_{01}
# plot_h1(1, 0, 0, 0, 0, 0)

# # T_{02} = P_{02}
# plot_h1(0, 1, 0, 0, 0, 0)

# # T_{03} = P_{03}
# plot_h1(0, 0, 1, 0, 0, 0)



### Example 1: Chebyshev Polynomial of order 1 with family 1
### Goal: Find the minimal sup-norm of P_{11} + a * P_{01}
### From a theoretical result, we know that a = -1/12 

### Plotting Chebyshev Polynomial of order 1 with family 1
# plot_h1(-1/12, 0, 0, 1, 0, 0, 7)

### Verifying the extremal points of the Chebyshev Polynomial
# extremal_val_h1(-1/12, 0, 0, 1, 0, 0, num_points=10, flag=3, level=7)


### Example 2: Chebyshev Polynomial of order 1 with family 2
### Goal: Find the minimal sup-norm of P_{12} + a * P_{02} 
### We find the optimal value through a grid-search procedure.

# # This indicates we evaluate 20 points each time.
# num_points = 20

### The various sets of start/end values that we used for finding the 
### approximate best coefficients by shrinking the range.
### Only one of those should be active.

# # Initial start and end values
# start = float(-2 * beta(1) / beta(0))
# end = 0
# level = 7

# # Second start and end values
# start = -0.0748538
# end = -0.0561403
# level = 7

# # Third start and end values
# start = -0.0630347
# end = -0.0610649
# level = 7

# # Fourth start and end values
# start = -0.0621016
# end = -0.0618943
# level = 10

# # Fifth start and end values
# start = -0.0619488
# end = -0.0619270
# level = 10

# # Sixth start and end values
# start = -0.0619350
# end = -0.0619327
# level = 10

# # Seventh start and end values
# start = -0.0619340
# end = -0.0619337
# level = 10

### Here we perform the grid search.

# # Print the maximal norm at each coefficient
# max_h1_val_family(start, end, num_points, 2, level)

### After running the above experiments, we find that the optimal 
### coefficient is around -0.0619339

# best_coeff = -0.0619339

### Well, something we noted in the process is that the sup-norm seems
### to be achieved at the point where the function has the same 
### norm at F_0F_1q_2 and at q_1/q_2. 
### Based on this conjecture, we can solve a linear equation that gives
### the exact value of a.

# addr = [0, 1, 2]
# best_coeff = -(p_jk_addr_list(addr, 1, 2) + beta(1))  \
#     / (p_jk_addr_list(addr, 0, 2) + beta(0))

### Plotting Chebyshev Polynomial of order 1 with family 2
# plot_h1(0, best_coeff, 0, 0, 1, 0, 7)

### Verifying the extremal points of the Chebyshev Polynomial
# extremal_val_h1(0, best_coeff, 0, 0, 1, 0, \
#     num_points=10, flag=3, level=7)



### Example 3: Chebyshev Polynomial of order 1 with family 3
### Goal: Find the minimal sup-norm of P_{13} + a * P_{03} 
### We find the optimal value through a grid-search procedure.

# # This indicates we evaluate 20 points one time.
# num_points = 20

### The various sets of start/end values that we used for finding the 
### approximate best coefficients by shrinking the range.
### Only one of those should be active.

# # Initial start and end values
# start = float(-2 * alpha(2) / alpha(1))
# end = 0
# level = 7

# # Second set of start and end values
# start = -0.0315789
# end = -0.0245614
# level = 7

# # Third set of start and end values
# start = -0.0278855
# end = -0.0271468
# level = 7

# # Fourth set of start and end values
# # From here we start using 10 levels
# start = -0.0275744 
# end = -0.0274578
# level = 10

# # Fifth set of start and end values
# start = -0.0275068 
# end = -0.0274946
# level = 10

# # Fifth set of start and end values
# start = -0.0275035
# end = -0.0275023
# level = 10

# # Sixth set of start and end values
# start = -0.0275024
# end = -0.0275023
# level = 10

### Here we perform the grid search.

# # Print the maximal norm at each coefficient
# max_h1_val_family(start, end, num_points, 3, level)

### After running the above experiments, we find that the optimal 
### coefficient is around -0.02750235

# best_coeff = -0.02750235

### Plotting Chebyshev Polynomial of order 1 with family 3
# plot_h1(0, 0, best_coeff, 0, 0, 1, 7)

### Verifying the extremal points of the Chebyshev Polynomial
# extremal_val_h1(0, 0, best_coeff, 0, 0, 1, \
#     num_points=10, flag=3, level=7)


### Experiment 4: Chebyshev Polynomial of order 2 with family 1
### Goal: Find the minimal sup-norm of P_{21} + a * P_{11} + b * P_{01}
### We find the optimal value through a grid-search procedure.
 
# # This indicates we evaluate 500*500 points one time.
# num0 = 500
# num1 = 500

# # Initial start and end values
# start0 = -1/180
# end0 = 1/180
# start1 = -1/10
# end1 = 1/30
# level = 7

### Here we perform the grid search.

# max_h2_val_family(start0, end0, num0, start1, end1, num1, 10, 1, level)

### After running the above experiments, we find that the optimal 
### coefficient is around a = -0.033199732798931195, 
### b = 0.0008795368514807394

### Well, something we noted in the process is that the sup-norm seems
### to be achieved at the point where the function has the same 
### norm at q_0, q_1, q_2 (max) and F1q2 (min)
### Based on this conjecture, we can solve a linear equation that gives
### the exact value of a and b

# addr = [1, 2]
# print(p_jk_addr_list(addr, 2, 1))
# print(p_jk_addr_list(addr, 1, 1))

# best_coeff_0 = 1 / 1125
# best_coeff_1 = - 1 / 30

### Here we generate the polynomialo with the best coefficients

# best_coeff_arr = np.zeros((3, 3))
# best_coeff_arr[0, 0] = best_coeff_0
# best_coeff_arr[1, 0] = best_coeff_1
# best_coeff_arr[2, 0] = 1


### Plotting Chebyshev Polynomial of order 2 with family 1
# plot_hn(best_coeff_arr, 7)

### Verifying the extremal points of the Chebyshev Polynomial
# extremal_val_hn(best_coeff_arr, 10, flag=3, level=7)


### Experiment 5: Chebyshev Polynomial of order 2 with family 2
### Goal: Find the minimal sup-norm of P_{22} + a * P_{12} + b * P_{02}
### We find the optimal value through a grid-search procedure.
 
# This indicates we evaluate 1000*1000 points one time.
# num0 = 200
# num1 = 200

# Initial start and end values
# start0 = -398/91125
# end0 = 751/182250
# start1 = -1532/2025
# end1 = 1568/2025
# level = 7

# Second set of start and end values
# start0 = 0.0004925638944843335 - 0.0001
# end0 = 0.0004925638944843335 + 0.0001
# start1 = -0.027122431073048348 - 0.002
# end1 = -0.027122431073048348 + 0.002
# level = 7

# Third set of start and end values
# start0 = 0.000467038368958808 - 0.0001
# end0 = 0.000467038368958808 + 0.0001
# start1 = -0.026215524166141438 - 0.002
# end1 = -0.026215524166141438 + 0.002
# level = 7

## Here we perform the grid search.

# max_h2_val_family(start0, end0, num0, start1, end1, num1, 10, 2, level)

### After running the above experiments, we find that the optimal 
### coefficient is around a = -0.033199732798931195, 
### b = 0.0008795368514807394 

### Well, something we noted in the process is that the sup-norm seems
### to be achieved at the point where the function has the same 
### norm at q_0, q_1, q_2 (max) and F1q2 (min)
### Based on this conjecture, we can solve a linear equation that gives
### the exact value of a and b

# addr = [1, 2]
# print(p_jk_addr_list(addr, 2, 1))
# print(p_jk_addr_list(addr, 1, 1))

# FIRST SET OF BEST VALUES
# best_coeff_0 = 0.0004925638944843335
# best_coeff_1 = -0.027122431073048348
# best_result = 8.185607160469618e-05

# best_coeff_0 = 0.000467038368958808
# best_coeff_1 = -0.026215524166141438
# best_result_2 = 7.675450852596407e-05

### Here we generate the polynomial with the best coefficients

# best_coeff_arr = np.zeros((3, 3))
# best_coeff_arr[0, 1] = best_coeff_0
# best_coeff_arr[1, 1] = best_coeff_1
# best_coeff_arr[2, 1] = 1


### Plotting Chebyshev Polynomial of order 2 with family 1
# plot_hn(best_coeff_arr, 7)

### Verifying the extremal points of the Chebyshev Polynomial
# extremal_val_hn(best_coeff_arr, 10, flag=3, level=7)


### Experiment 6: Chebyshev Polynomial of order 2 with family 3
### Goal: Find the minimal sup-norm of P_{23} + a * P_{13} + b * P_{03}
### We find the optimal value through a grid-search procedure.
 
# This indicates we evaluate 500*500 points one time.
# num0 = 1000
# num1 = 1000

# Initial start and end values
# start0 = -7/3000
# end0 = 67/27000
# start1 = -29/300
# end1 = 7/100
# level = 7

# Second set of start and end values
# start0 = 0.0001535980424869311 - 0.00002
# end0 = 0.0001535980424869311 + 0.00002
# start1 = -0.01541875208541875 - 0.002
# end1 = -0.01541875208541875 + 0.002

# level = 7

## Here we perform the grid search.

# max_h2_val_family(start0, end0, num0, start1, end1, num1, 10, 3, level)

### After running the above experiments, we find that the optimal 
### coefficient is around a = -0.033199732798931195, 
### b = 0.0008795368514807394 

### Well, something we noted in the process is that the sup-norm seems
### to be achieved at the point where the function has the same 
### norm at q_0, q_1, q_2 (max) and F1q2 (min)
### Based on this conjecture, we can solve a linear equation that gives
### the exact value of a and b

# addr = [1, 2]
# print(p_jk_addr_list(addr, 2, 1))
# print(p_jk_addr_list(addr, 1, 1))

# FIRST SET OF BEST VALUES
# best_coeff_0 = 0.0001535980424869311
# best_coeff_1 = -0.01541875208541875
# best_result = 7.3784509667491625e-06

# best_coeff_0 = 0.00015121566010454872
# best_coeff_1 = -0.015236569903236567

### Here we generate the polynomial with the best coefficients

# best_coeff_arr = np.zeros((3, 3))
# best_coeff_arr[0, 2] = best_coeff_0
# best_coeff_arr[1, 2] = best_coeff_1
# best_coeff_arr[2, 2] = 1


### Plotting Chebyshev Polynomial of order 2 with family 1
# plot_hn(best_coeff_arr, 7)

### Verifying the extremal points of the Chebyshev Polynomial
# extremal_val_hn(best_coeff_arr, 10, flag=3, level=7)


###---------------  VI. Full Chebyshev Polynomials -----------------###


### Experiment 1: Full Chebyshev Polynomial of order 1 with family 1
### Goal: Find the minimal sup-norm of P_{11} + a * P_{01} + b * P_{02} 
###     + c * P_{03}
### We find the optimal value through a grid-search procedure.

# # This indicates we evaluate 500*500 points one time.
# num1 = 100
# num2 = 100
# num3 = 100

### Initial start and end values
# start1 = -1/6
# end1 = 1/6
# start2 = -1/3
# end2 = 1
# start3 = -1/3
# end3 = 1/3
# level = 7

### First best 0.03535353535353536


### Second start and end values
# start1 = 0.03535353535353536 - 0.01
# end1 = 0.03535353535353536 + 0.01
# start2 = 0.34006734006734 + 0.05
# end2 = 0.34006734006734 - 0.05
# start3 = -0.0033670033670033517 - 0.01
# end3 = -0.0033670033670033517 + 0.01
# level = 7

### Second best 0.03343434343434344


### Third start and end values
# start1 = 0.03343434343434344 - 0.0005
# end1 = 0.03343434343434344 + 0.0005
# start2 = 0.33350168350168347 + 0.001
# end2 = 0.33350168350168347 - 0.001
# start3 = -3.367003367001943e-05 - 0.0002
# end3 = -3.367003367001943e-05 + 0.0002
# level = 7

### Here we perform the grid search.

# max_h1_full_family(start1, end1, num1, start2, end2, num2, start3, \
#     end3, num3, 10, 1, level)

### THIRD BEST COEFFICIENTS
# best_coeff_1 = 0.03332828282828283
# best_coeff_2 = 0.3333299663299663
# best_coeff_3 = 6.734006734149093e-07

### Well, we note that the above values are very close to a = 1/30, 
### b=1/3, c=0. We can verify this does give the desired properties.

# best_coeff_1 = 1/30
# best_coeff_2 = 1/3
# best_coeff_3 = 0


### Plotting Full Chebyshev Polynomial of order 1 with family 1
# plot_h1(best_coeff_1, best_coeff_2, best_coeff_3, 1, 0, 0, 7)

### Verifying the extremal points of the Chebyshev Polynomial
# extremal_val_h1(best_coeff_1, best_coeff_2, best_coeff_3, 1, 0, 0, 10, 3, 7)



### Experiment 2: Full Chebyshev Polynomial of order 1 with family 2
### Goal: Find the minimal sup-norm of P_{12} + a * P_{01} + b * P_{02} 
###     + c * P_{03}
### We find the optimal value through a grid-search procedure.

# # This indicates we evaluate 500*500 points one time.
num1 = 100 + 1
num2 = 100 + 1
num3 = 100 + 1

### Initial start and end values
# start1 = -2/45
# end1 = 2/45
# start2 = -12/45
# end2 = 4/45
# start3 = -4/45
# end3 = 4/45
# level = 7

### First best 0.012444444444444439

### Second start and end values
# start1 = -0.011555555555555555 - 0.005
# end1 = -0.011555555555555555 + 0.005
# start2 = -0.08888888888888888 + 0.01
# end2 = -0.08888888888888888 - 0.01
# start3 = 0
# end3 = 0
# level = 7

### Second best 0.012044444444444437


### Third start and end values
# start1 = -0.011955555555555556 - 0.0002
# end1 = -0.011955555555555556 + 0.0002
# start2 = -0.08888888888888888
# end2 = -0.08888888888888888
# start3 = 0
# end3 = 0
# level = 7

### Third best 0.012000444444444439


### Here we perform the grid search.

# max_h1_full_family(start1, end1, num1, start2, end2, num2, start3, \
#     end3, num3, 10, 2, level)

### THIRD BEST COEFFICIENTS
# best_coeff_1 = -0.011999555555555557
# best_coeff_2 = -0.08888888888888888
# best_coeff_3 = 0

### Well, we note that the above values are very close to a = -3/250, 
### b=-4/45, c=0. We can verify this does give the desired properties.

# best_coeff_1 = -3/250
# best_coeff_2 = -4/45
# best_coeff_3 = 0


### Plotting Full Chebyshev Polynomial of order 1 with family 1
# plot_h1(best_coeff_1, best_coeff_2, best_coeff_3, 0, 1, 0, 7)

### Verifying the extremal points of the Chebyshev Polynomial
# extremal_val_h1(best_coeff_1, best_coeff_2, best_coeff_3, 0, 1, 0, 10, 3, 7)

### Some preliminary computations for finding the best values
# addr = [1, 2]
# print(p_jk_addr_list(addr, 1, 2))
# print(p_jk_addr_list(addr, 0, 1))
# print(p_jk_addr_list(addr, 0, 2))
# print(p_jk_addr_list(addr, 0, 3))


### Experiment 3: Full Chebyshev Polynomial of order 1 with family 3
### Goal: Find the minimal sup-norm of P_{13} + a * P_{01} + b * P_{02} 
###     + c * P_{03}
### We find the optimal value through a grid-search procedure.

# # This indicates we evaluate 500*500 points one time.
num1 = 1
num2 = 1
num3 = 100 + 1

### Initial start and end values
# start1 = -1/60
# end1 = 1/60
# start2 = -1/15
# end2 = 1/15
# start3 = -1/15
# end3 = 0
# level = 7

### First best 0.002999999999999999

### Second start and end values
# start1 = 0 - 0.002
# end1 = 0 + 0.002
# start2 = 0 - 0.01
# end2 = 0 + 0.01
# start3 = -0.027333333333333334 - 0.002
# end3 = -0.027333333333333334 + 0.002
# level = 7

### Second best 0.002919999999999999


### Third start and end values
# start1 = 0 - 0.0002
# end1 = 0 + 0.0002
# start2 = 0 - 0.0004
# end2 = 0 + 0.0004
# start3 = -0.027493333333333335 - 0.0002
# end3 = -0.027493333333333335 + 0.0002
# level = 7

### Third best 0.002916007389866667

### Fourth start and end values
# start1 = 0
# end1 = 0
# start2 = 0
# end2 = 0
# start3 = -0.027503333333333334 - 0.0001
# end3 = -0.027503333333333334 + 0.0001
# level = 7

### Fourth best 0.002915651037866667

### Fifth start and end values
# start1 = 0
# end1 = 0
# start2 = 0
# end2 = 0
# start3 = -0.027501333333333336 - 0.00002
# end3 = -0.027501333333333336 + 0.00002
# level = 7

### Fifth best 0.0029155084970666667

### Sixth start and end values
# start1 = 0
# end1 = 0
# start2 = 0
# end2 = 0
# start3 = -0.027502533333333336 - 0.0000002
# end3 = -0.027502533333333336 + 0.0000002
# level = 7

### Sixth best 0.002915479999999998

### Seventh start and end values
start1 = 0
end1 = 0
start2 = 0
end2 = 0
start3 = -0.027502373333333337 - 0.00000001
end3 = -0.027502373333333337 + 0.00000001
level = 10

### Here we perform the grid search.

max_h1_full_family(start1, end1, num1, start2, end2, num2, start3, \
    end3, num3, 10, 3, level)

### THIRD BEST COEFFICIENTS
best_coeff_1 = 0
best_coeff_2 = 0
best_coeff_3 = -0.027502373333333337

### Well, we note that the above values are the same as the special 
# chebyshev polynomial of order 1 of family 3, so we share the same
# graphs.


### Plotting Full Chebyshev Polynomial of order 1 with family 1
# plot_h1(best_coeff_1, best_coeff_2, best_coeff_3, 0, 0, 1, 7)

# ### Verifying the extremal points of the Chebyshev Polynomial
# extremal_val_h1(best_coeff_1, best_coeff_2, best_coeff_3, 0, 0, 1, 10, 3, 7)


###-------------- V. Finding Optimality Conditions ---------------###

### Experiment 1: P_{03} / P_{11}

# # Initiaize the parameters
# num_points = 10

# # Generate the values of the polynomial
# T = generate_T(7, 2, False)
# P01_seq = T[0, :, 0]
# P02_seq = T[1, :, 0]
# P03_seq = T[2, :, 0]
# P11_seq = T[0, :, 1]
# P12_seq = T[1, :, 1]
# P13_seq = T[2, :, 1]

# # Evaluate the target polynomial
# div_P03_P11 = np.divide(P03_seq, P11_seq)

# # Find the extremal value
# extremal_val(div_P03_P11, num_points, flag=1, level=7)

# # Plot the resulting polynomial
# plot_general(div_P03_P11, level=7)


### Experiment 2: f_{10} + f_{11} + f_{12}

# # Initiaize the parameters
# num_points = 10

# # Generate the values of the polynomial
# W = generate_W(7, 2, False)
# f00_seq = W[0, :, 0]
# f01_seq = W[1, :, 0]
# f02_seq = W[2, :, 0]
# f10_seq = W[0, :, 1]
# f11_seq = W[1, :, 1]
# f12_seq = W[2, :, 1]

# # Evaluate the target polynomial
# sum_f10_f11_f12 = f10_seq + f11_seq + f12_seq

# # Find the extremal value
# extremal_val(sum_f10_f11_f12, num_points, flag=1, level=7)

# # Plot the resulting polynomial
# plot_general(sum_f10_f11_f12, level=7)


### Experiment 3: P_{21} / P_{01}, P_{11} / P_{01}

# # Initiaize the parameters
# num_points = 10

# # Generate the values of the polynomial
# T = generate_T(7, 3, False)
# P01_seq = T[0, :, 0]
# P11_seq = T[0, :, 1]
# P21_seq = T[0, :, 2]

# # Evaluate the target polynomials
# div_P21_P01 = np.divide(P21_seq, P01_seq)
# div_P11_P01 = np.divide(P11_seq, P01_seq)

# # Find the extremal value
# extremal_val(div_P21_P01, num_points, flag=1, level=7)
# extremal_val(div_P11_P01, num_points, flag=1, level=7)

# # Plot the resulting polynomial
# plot_general(div_P21_P01, level=7)
# plot_general(div_P11_P01, level=7)


### Experiment 4: P_{23} / P_{03}, P_{13} / P_{03}

# # Initiaize the parameters
# num_points = 30

# # Generate the values of the polynomial
# T = generate_T(7, 3, False)
# P03_seq = T[2, :, 0]
# P13_seq = T[2, :, 1]
# P23_seq = T[2, :, 2]

# # Evaluate the target polynomials
# div_P23_P03 = np.divide(P23_seq, P03_seq)
# div_P13_P03 = np.divide(P13_seq, P03_seq)

# # Find the extremal value
# # extremal_val(div_P23_P03, num_points, flag=1, level=7)
# extremal_val(div_P13_P03, num_points, flag=1, level=7)

# # Plot the resulting polynomial
# plot_general(div_P23_P03, level=7)
# plot_general(div_P13_P03, level=7)

# print(p_jk_addr_list([0, 1], 2, 3))
# print(p_jk_addr_list([0, 1], 1, 3))
# print(p_jk_addr_list([0, 1], 0, 3))
# print(p_jk_addr_list([1, 1], 2, 3))
# print(p_jk_addr_list([1, 1], 1, 3))
# print(p_jk_addr_list([1, 1], 0, 3))

# print(gamma(0))
# print(gamma(1))
# print(beta(2))
# print(p_jk_addr_list([0, 1], 1, 1))
# print(p_jk_addr_list([0, 1], 0, 1))
# print(p_jk_addr_list([0, 1], 0, 2))
# print(p_jk_addr_list([0, 1], 0, 3))
# print(p_jk_addr_list([0, 1], 1, 2))
# print(p_jk_addr_list([0, 1], 0, 2))



