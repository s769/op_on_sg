import sys

import numpy as np
import gmpy2 as gm

from util import address_from_index
from recursions import alpha, beta
from monomials import p_jk, generate_W, generate_T, generate_norm_W, generate_norm_T, norm_p_jk, f_lkFiqn, norm_f_jk
from plotting import plot_monomial, plot_easy_basis, plot_op, plot_general
from chebyshev import eval_poly_h1, plot_h1, extremal_val, extremal_val_h1, max_h1_val_family
np.set_printoptions(threshold=sys.maxsize)

# An outline of the functionalities are listed here.
# Section I. Plotting
# 1. Plotting the easy basis/monomial/SOP
# 2. Plotting a general polynomial


#### Section I. Plotting

######## 1. Code for Plotting the Easy Basis/Monomial/SOP ########

### Parameters
# num represent the highest power to be plotted
# k represent the type of monomial {0,1,2} -> {1, 2, 3}
# num = 3
# k = 0


# 1-a. To plot the monomials, use plot_monomial
#plot_monomial(num, k)


# 1-b. To plot the easy basis, use plot_easy_basis
#plot_easy_basis(num, k)


# 1-c. This plots the Sobolov Orthogonal Polynomials
#plot_SOP(num, k)


######## 2. Code for Plotting a general polynomial in H1 ########

# Example: Compute P_{13} + 1/2*P_{03}
#plot_h1(0, 0, 1/2, 0, 0, 1, 7)



### Section II. Maximal, Minimal, Sup-norms

######## 1. Code for finding sup-norm of general polynomials ########

### Example A: Sup-norm of f_{j0}
# We can see from the following code that f_{j0} always take the maximum
# norm at point F_0F_1q_2=F_0F_2q_1

# # Initiaize the parameters
# j = 3
# num_points = 10

# # Generate the values of the polynomial
# W = generate_W(7, j+1, False)
# fj1_val = W[0, :, j]

# # Find the extremal value
# extremal_val(fj1_val, num_points, flag=3, level=7)


### Example B: Sup-norm of f_{j+1, k} / f_{j, k} 
# We can see that the sup-norm is achieved at F_{1}q_{2}

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


### Example C: Sup norm of P_{03} / P_{11}
# The result is not very clear - seems that a lot of points have
# similar values.

# # Initiaize the parameters
# j = 0
# num_points = 10

# # Generate the values of the polynomial
# T = generate_T(7, j+2, False)
# P03_seq = T[2, :, 0]
# P11_seq = T[0, :, 1]
# prop = np.abs(np.divide(P03_seq, P11_seq))

# # Find the extremal value
# extremal_val(prop, num_points, flag=3, level=7)


### Actual work 1: Sup norm of P_{13} + a * P_{03}

# Initialization of the coefficients

### Here we find the optimal value theore
# # This indicates we evaluate 20 points one time.
# num_points = 20

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

# # Print the maximal norm at each coefficient
# max_h1_val_family(start, end, num_points, 3, level)

# # This implies the optimal coefficient is around -0.02750235
# best_coeff = -0.02750235

# # We finally print the graph of the final polynomial we obtain
# plot_h1(0, 0, best_coeff, 0, 0, 1, 7)


### Actual work 2: Sup norm of P_{12} + a * P_{02}

### Here we find the optimal value theore
# # This indicates we evaluate 20 points one time.
# num_points = 20

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

# # Sixth start and end values
# start = -0.0619340
# end = -0.0619337
# level = 10

# # Print the maximal norm at each coefficient
# max_h1_val_family(start, end, num_points, 2, level)

# # This implies the optimal coefficient is around -0.0619339
# best_coeff = -0.0619339

# Well, something we noted in the process is that the sup-norm seems
# to be achieved at the point where the function has the same 
# norm at F_0F_1q_2 and at q_1/q_2. 

# Based on this conjecture, we can solve a linear equation that gives
# the exact value of a:
addr = [0, 1, 2]
addr = np.flip(addr)
addr = ''.join(str(int(x)) for x in addr)
best_coeff = -(p_jk(addr, 1, 2) + beta(1)) / (p_jk(addr, 0, 2) + beta(0))
print(float(best_coeff))
print(float(beta(1) + best_coeff * beta(0)))

# # We finally print the graph of the final polynomial we obtain
# plot_h1(0, best_coeff, 0, 0, 1, 0, 7)

#extremal_val_h1(10, 0, best_coeff, 0, 0, 1, 0, 3, 7)

# ###
# addr = [0, 1, 2]
# addr = np.flip(addr)
# addr = ''.join(str(int(x)) for x in addr)

# print(p_jk(addr, 1, 2))
# print(p_jk(addr, 0, 2))
# # a = -(p_jk(addr, 1, 2) - beta(1)) / (p_jk(addr, 1, 2) - beta(0))
# # or negative
# a = -(p_jk(addr, 1, 2) + beta(1)) / (p_jk(addr, 0, 2) + beta(0))
# print(a)