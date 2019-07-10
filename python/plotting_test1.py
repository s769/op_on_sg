import sys
import numpy as np

from monomials import generate_W
from plotting import plot_monomial, plot_easy_basis, plot_SOP
from util import address_from_index

np.set_printoptions(threshold=sys.maxsize)

### I. Code for Plotting the Easy Basis/Monomial/SOP ###

# 1. To plot the monomials, use plot_monomial
# num represent the highest power to be plotted
# k represent the type of monomial {0,1,2} -> {1, 2, 3}
num = 8
k = 0

#plot_monomial(num, k)

# 2. To plot the easy basis, use plot_easy_basis
#plot_easy_basis(num, k)

# 3. This plots the Sobolov Orthogonal Polynomials
# num represent the highest power to be plotted
# k represent the type of monomial {0,1,2} -> {1, 2, 3}
#plot_SOP(num, k)



### II. Verifying the Sup norm of the easy basis f_{jk} ###

# Since the f_{jk}'s are just rotations of each other,
# we can just assume k = 0.

# Order of the easy basis j
j = 5

# Fetch the values for f_{jk}
W = generate_W(7, j+1, False)
f10_seq = W[0, :, j]

# Sort the values on all the points 
if (j % 2 == 0):
    ind_arrange = np.flip(np.argsort(f10_seq))
else:
    ind_arrange = np.argsort(f10_seq)

# Find the address of the maximal points
for i in range(10):
    temp_add = address_from_index(7, ind_arrange[i]+1)
    print("Address ", i)
    print(temp_add)
    print("Index", i)
    print(ind_arrange[i])
    print("Value ", i)
    print(f10_seq[ind_arrange[i]])

# We can see from the printed statements that the sup norm
# is always achieved at F_0F_1q_2=F_0F_2q_1







