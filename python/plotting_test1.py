import sys

import numpy as np
import gmpy2 as gm

from util import address_from_index
from recursions import alpha
from monomials import generate_W, generate_T, generate_norm_W, generate_norm_T, norm_p_jk, f_lkFiqn, norm_f_jk
from plotting import plot_monomial, plot_easy_basis, plot_op, plot_general
from chebyshev import eval_poly_h1, plot_h1_family, plot_h1, max_h1_family, max_h1_val, max_h1_val_family

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



######## 2. Code for Plotting a general polynomial ########

# Example: Compute P_{13} + a*P_{03}

# Hyperparameters
#a = 1 / 2

# Maximum order of the basis we're working with
#j = 1

# # Fetch the values for P
# T = generate_T(7, j+2, False)
# P13_seq = T[2, :, 1]
# P03_seq = T[2, :, 0]

# # print()
# # print()
# # print()

#plot_general(prop, level)

######## II. Verifying the Sup norm of the easy basis f_{jk} ########

# Since the f_{jk}'s are just rotations of each other,
# we can just assume k = 0.

# # Order of the easy basis j
# j = 5

# # Fetch the values for f_{jk}
# W = generate_W(7, j+1, False)
# fj0_seq = W[0, :, j]

# # Sort the values on all the points 
# if (j % 2 == 0):
#     ind_arrange = np.flip(np.argsort(f10_seq))
# else:
#     ind_arrange = np.argsort(f10_seq)

# # Find the address of the maximal points
# for i in range(10):
#     temp_add = address_from_index(7, ind_arrange[i]+1)
#     print("Address ", i)
#     print(temp_add)
#     print("Index", i)
#     print(ind_arrange[i])
#     print("Value ", i)
#     print(fj0_seq[ind_arrange[i]])

# We can see from the printed statements that the sup norm
# is always achieved at F_0F_1q_2=F_0F_2q_1




######## 3. Verifying the sup norm of f_{j+1,k} / f_{j, k} ########

# Again we can assume k = 0.

# # Order of the easy basis j (>=0)
# j = 7

# # Fetch the values for f_{jk}
# W = generate_W(7, j+2, False)
# fj0_seq = W[0, :, j]
# fj1_0_seq = W[0, :, j+1]

# # print()
# # print()
# # print()
# # print(fj0_seq.size)
# # print(fj1_0_seq.size)

# prop = np.abs(np.divide(fj1_0_seq, fj0_seq))

# ind_arrange_1 = np.flip(np.argsort(prop))
# ind_arrange_2 = np.argsort(prop)

# for i in range(10):
#     temp_add = address_from_index(7, ind_arrange_1[i]+1)
#     print("Address ", i)
#     print(temp_add)
#     print("Index", i)
#     print(ind_arrange_1[i])
#     print("Value ", i)
#     print(prop[ind_arrange_1[i]])

# for i in range(10):
#     temp_add = address_from_index(7, ind_arrange_2[i]+1)
#     print("Address ", (i+10))
#     print(temp_add)
#     print("Index", (i+10))
#     print(ind_arrange_2[i])
#     print("Value ", (i+10))
#     print(prop[ind_arrange_2[i]])


######## 4. Verifying the value of \partial_{n}p_{jk} ########

# Again we can assume k = 0.

# Order of the easy basis j (>=0)
# j = 4

# Fetch the values for \partial_{n}p_{jk}
# norm_W = generate_norm_W(1, j+2, False)
# norm_f00 = norm_W[0, :, j]
# print(norm_f00)

# Fetch the values for \partial_{n}p_{jk}
# norm_T = generate_norm_T(1, j+2, False)
# norm_P01 = norm_T[2, :, j]
# print()
# print()
# print(norm_P01)


### 5. Verifying the sup norm of P_{03} / P_{11} ###

#Again we can assume k = 0.

#level = 7

# Fetch the values for P_{03} and P_{11}
#T = generate_T(level, 2, False)
#P03_seq = T[2, :, 0]
#P11_seq = T[0, :, 1]

# print()
# print()
# print()
# print(fj0_seq.size)
# print(fj1_0_seq.size)

#prop = np.divide(P03_seq, P11_seq)

#plot_general(prop, level)



# ind_arrange_1 = np.flip(np.argsort(prop))
# ind_arrange_2 = np.argsort(prop)

# for i in range(10):
#     temp_add = address_from_index(level, ind_arrange_1[i]+1)
#     print("Address ", i)
#     print(temp_add)
#     print("Index", i)
#     print(ind_arrange_1[i])
#     print("Value ", i)
#     print(prop[ind_arrange_1[i]])

# for i in range(10):
#     temp_add = address_from_index(level, ind_arrange_2[i]+1)
#     print("Address ", (i+10))
#     print(temp_add)
#     print("Index", (i+10))
#     print(ind_arrange_2[i])
#     print("Value ", (i+10))
#     print(prop[ind_arrange_2[i]])


### 3. Verifying the sup norm of f_{j+1,k} / f_{j, k} ###

# Again we can assume k = 0.

# # Order of the easy basis j (>=0)
# j = 7

# # Fetch the values for f_{jk}
# T = generate_T(1, j+2, False)
# P11_seq = T[0, :, 6]

# print()
# print()
# print(P11_seq)


#plot_h1(0, 0, float(-alpha(2) / (2*alpha(1))), 0, 0, 1, 7)
#start = float(-2*alpha(2)/alpha(1))
#end = 0

#start = float(-alpha(2)/alpha(1) - alpha(2)/(4*alpha(1)))
#end = float(-alpha(2)/alpha(1) + alpha(2)/(4*alpha(1)))

start = -0.02757
end = -0.02747


#start = -0.0275068
#end = -0.0275015

#plot_h1_family(start, end, 10, 3, 7)
# T = eval_poly_h1(0, 0, -0.027, 0, 0, 1, 7)
# T = np.abs(T)
# ind_arrange_1 = np.flip(np.argsort(T))

# for i in range(10):
#     temp_add = address_from_index(7, ind_arrange_1[i]+1)
#     print("Address ", i)
#     print(temp_add)
#     print("Index", i)
#     print(ind_arrange_1[i])
#     print("Value ", i)
#     print(T[ind_arrange_1[i]])

#max_h1_val(-0.0275068, 3, 10, 9)

#new_range = np.linspace(start, end, 10)
#max_h1_family(start, end, 20, 3, 7)
#max_h1_val_family(start, end, 20, 3, 10)

# addr = address_from_index(7, 3 ** 8)
# print(addr)
# addr = np.flip(addr)
# addr = ''.join(str(int(x)) for x in addr)

# print(norm_p_jk(addr, 1, 3) - 0.0275013 * norm_p_jk(addr, 0, 3))
# print(norm_p_jk(addr, 1, 3))
# print(norm_p_jk(addr, 0, 3))

# #max_h1_val(-0.0275013, 3, 10, 7)

plot_h1(-1/12,0,0, 1, 0, 0, 7)
#plot_h1(0,0,-0.0275013, 0, 0, 1, 7)

#print(alpha(3))