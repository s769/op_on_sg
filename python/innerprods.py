import numpy as np
from recursions import alpha, beta, gamma, eta, ap, alpha_array,\
    beta_array, gamma_array, eta_array, ap_array
import functools
# import sympy as sp
# from sympy import Rational as Rat
import gmpy2 as gm
'''
This file contains functions used to compute the L2 inner product 
between the monomial baisis P_jk. The algorithms for computing the inner 
products are derived from the Kasso, Tuley paper OP on SG.

There are some slight corrections to the algorithms mentioned in the 
paper that are implemented here.

'''

# def inner_j1k1(j, k, energy=False):
#   if energy:
#     res = 2*alpha(k)*eta(j)
#     ms = min(j, k)
#     s1 = 0
#     for l in range(j-ms, j+1):
#       s1 += alpha(j-l)*eta(k+l+1) - alpha(k+l+1)*eta(j-l)
#     s2 = 0
#     for l in range(j-ms-1, j):
#       s2 += alpha(j-l-1)*eta(k+l+1) - alpha(k+l+1)*eta(j-l-1)

#     return res + 2*(s1-s2)


#   ms = min(j, k)
#   s1 = 0
#   for l in range(j-ms, j+1):
#     s1 += alpha(j-l)*eta(k+l+1) - alpha(k+l+1)*eta(j-l)
#   s2 = 0
#   for l in range(j-ms-1, j):
#     s2 += alpha(j-l-1)*n(k+l) - alpha(k+l)*n(j-l-1)

#   return 2*(s1+s2)


def inner0_j1k1(j, k):
    '''
    Calculates the L2 inner product <P_j1, P_k1>

    Args:
        j, k: indices for the monomials P_j1, P_k1

    Returns:
        L2 inner product <P_j1, P_k1>

    '''

    ms = min(j, k)
    s1 = 0
    for l in range(j-ms, j+1):
        s1 += alpha(j-l)*eta(k+l+1) - alpha(k+l+1)*eta(j-l)
    return 2*s1

# def inner_j2k2(j, k, energy=False):
#   if energy:
#     res = -2*beta(k)*alpha(j)
#     s1 = 0
#     ms = min(j, k)
#     for l in range(j-ms, j+1):
#       s1 += beta(j-l)*alpha(k+l+1) - beta(k+l+1)*alpha(j-l)
#     s2 = 0
#     for l in range(j-ms-1, j):
#       s2 += beta(j-l-1)*alpha(k+l+1) - beta(k+l+1)*alpha(j-l-1)

#     return res - 2*(s1-s2)

#   ms = min(j, k)
#   s1 = 0
#   for l in range(j-ms, j+1):
#     s1 += beta(j-l)*alpha(k+l+1) - beta(k+l+1)*alpha(j-l)
#   s2 = 0
#   for l in range(j-ms-1, j):
#     s2 += beta(j-l-1)*alpha(k+l) - beta(k+l)*alpha(j-l-1)

#   return -2*(s1+s2)


def inner0_j2k2(j, k):
    '''
    Calculates the L2 inner product <P_j2, P_k2>

    Args:
        j, k: indices for the monomials P_j2, P_k

    Returns:
        L2 inner product <P_j2, P_k2>
    '''
    ms = min(j, k)
    s1 = 0
    for l in range(j-ms, j+1):
        s1 += beta(j-l)*alpha(k+l+1) - beta(k+l+1)*ap(j-l)
    return -2*s1

# def inner_j3k3(j, k, energy=False):
#   if energy:
#     res = 6*gamma(k)*eta(k)
#     ms = min(j, k)
#     s1 = 0
#     for l in range(j-ms, j+1):
#       s1 += alpha(j-l+1)*eta(k+l+2) - alpha(k+l+2)*eta(j-l+1)
#     s2 = 0
#     for l in range(j-ms-1, j):
#       s2 += alpha(j-l)*eta(k+l+2) - alpha(k+l+2)*eta(j-l)

#     return res + 18*(s1-s2)

#   ms = min(j, k)
#   s1 = 0
#   for l in range(j-ms, j+1):
#     s1 += alpha(j-l+1)*eta(k+l+2) - alpha(k+l+2)*eta(j-l+1)
#   s2 = 0
#   for l in range(j-ms-1, j):
#     s2 += alpha(j-l)*eta(k+l+1) - alpha(k+l+1)*eta(j-l)

#   return 18*(s1+s2)


def inner0_j3k3(j, k):
    '''
    Calculates the L2 inner product <P_j3, P_k3>

    Args:
        j, k: indices for the monomials P_j3, P_k3

    Returns:
        L2 inner product <P_j3, P_k3>
    '''
    ms = min(j, k)
    s1 = 0
    for l in range(j-ms, j+1):
        s1 += alpha(j-l+1)*eta(k+l+2) - alpha(k+l+2)*eta(j-l+1)
    return 18*s1

# def inner_j1k2(j, k, energy=False):
#   if energy:
#     res = 2*beta(k)*eta(k)
#     s1 = 0
#     for l in range(j+1):
#       s1 += alpha(j-l)*alpha(k+1+l) + beta(k+1+l)*eta(j-l)
#     s2 = 0
#     for l in range(j):
#       s2 += alpha(j-l-1)*alpha(k+1+l) + beta(k+1+l)*eta(j-l-1)

#     return res - 2*(s1-s2)

#   s1 = 0
#   for l in range(j+1):
#     s1 += alpha(j-l)*alpha(k+l+1) + beta(k+l+1)*eta(j-l)
#   s2 = 0
#   for l in range(j):
#     s2 += alpha(j-l-1)*alpha(k+l) + beta(k+l)*eta(j-l-1)

#   return -2*(s1+s2)


def inner0_j1k2(j, k):
    '''
    Calculates the L2 inner product <P_j1, P_k2>

    Args:
        j, k: indices for the monomials P_j1, P_k2

    Returns:
        L2 inner product <P_j1, P_k2>

    '''
    s1 = 0
    for l in range(j+1):
        s1 += alpha(j-l)*alpha(k+l+1) + beta(k+l+1)*eta(j-l)
    return -2*s1


# This function is used to symmetrize an upper triangular matrix.
# This function is used when creating Gram Matrices for the inner products.ArithmeticError
def symmetrize(arr):
    #return arr + arr.T - arr.multiply_elementwise(eye(arr.rows))
    return arr + arr.T - np.diag(np.diag(arr))

# This function takes a list/array of integers and outputs the concatenation of the integers


def lis2str(lis):
    '''
    Convert a list of integers to an integer string using concatenation.

    Args:
        lis: list or np.array of integers

    Returns:
        a string which is the concatenation of the numbers in lis

    Example:
        lis2str(np.array([1, 2, 3]))
        >> '123'
        lis2str([012, 345, 678])
        >> '012345678'
    '''
    return ''.join(str(int(x)) for x in lis)


# The L2 inner products <P_j1, P_k3> and <P_j2, P_k3> are 0
def inner0_j1k3(j, k): return 0


def inner0_j2k3(j, k): return 0

'''
This is a dictionary mapping the values (i, i') to the L2 inner product 
    function for <P_ji, P_ki'>. This dictionary is used in the 
    construction of the Polynomial class.
'''
inner_dict = {(1, 1): inner0_j1k1, (2, 2): inner0_j2k2, (3, 3): inner0_j3k3,
              (1, 2): inner0_j1k2, (1, 3): inner0_j1k3, (2, 3): inner0_j2k3}


'''
vals_dict maps the values 1, 2, 3 to the functions alpha, beta, and 
    gamma. This dictionary, along with norm_dict are used in the computation 
    of the values and normal derivatives of the polynomials P_jk (k = 1, 2, 3) 
    on the boundary of SG. The functions are based on the Kasso, Tuley paper.

'''
vals_dict = {1: alpha, 2: beta, 3: gamma}


def dnpj2(j): return -alpha(j)


def dnpj3(j): return 3*eta(j+1)


norm_dict = {1: eta, 2: dnpj2, 3: dnpj3}
