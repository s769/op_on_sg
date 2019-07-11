import time
import copy

import numpy as np
import gmpy2 as gm
import tqdm

from recursions import mem, mem2, alpha, beta, gamma, eta, tau, zeros_gm
from util import address_from_index

'''
This file contains the functions that will compute the values of the 
    easy/monomial basis on a given level of SG. The recursive functions 
    are memoized (see recursions.py).
'''


@mem2
def big_recursion(j):
    '''
    This function computes the coefficients a_j, b_j, p_j, q_j found in 
        the Splines on Fractals paper.

    Args:
        j: index of coefficients

    Returns:
        4 x j+1 np.array of coefficients a_j, b_j, p_j, q_j from 0 to j
    '''
    # Initialize space for the coefficients
    p_arr = zeros_gm(j+1,1)
    q_arr = zeros_gm(j+1,1)
    a_arr = zeros_gm(j+1,1)
    b_arr = zeros_gm(j+1,1)

    # Initialize arrays for the 0-th term
    a_arr[0] = gm.mpq(7, 45)
    b_arr[0] = gm.mpq(4, 45)
    p_arr[0] = gm.mpq(2, 5)
    q_arr[0] = gm.mpq(1, 5)

    if j == 0:
        return np.vstack((a_arr, b_arr, p_arr, q_arr))
    
    # Main recursion
    for l in tqdm.tqdm(range(1, j+1)):
        # Implements equation (5.6) in Splines paper
        res3 = 0
        res4 = 0
        for k in range(l):
            p = p_arr[k]
            q = q_arr[k]
            a = a_arr[l-k-1]
            b = b_arr[l-k-1]
            res3 += (4*a + 3*b)*p + (a + 2*b)*q
            res4 += (2*a + 4*b)*p + (3*a + b)*q
        
        b = b_arr[l-1]
        p_arr[l] = -gm.mpq(2,5)*b - gm.mpq(1,5)*res3
        q_arr[l] = -gm.mpq(1,5)*b - gm.mpq(1,5)*res4

        # Implements equation (5.5) in Splines paper
        res1 = 0
        res2 = 0
        vec = zeros_gm(2,1)
        for k in range(l):
            p = p_arr[l-k]
            q = q_arr[l-k]
            a = a_arr[k]
            b = b_arr[k]
            res1 += (2*p + q)*(a+2*b)
            res2 += (p + 2*q)*(a + 2*b)
        
        vec[0] = gm.mpq(2,(3*(5**l - 1)))*res1
        vec[1] = gm.mpq(10,(3*(5**(l+1) - 1)))*res2
        coefs_inv = np.array([[gm.mpq(7,15), gm.mpq(-2,15)], [gm.mpq(4,15), gm.mpq(1,15)]])
        ab_arr = coefs_inv.dot(vec)
        a_arr[l] = ab_arr[0]
        b_arr[l] = ab_arr[1]
    
    return np.vstack((a_arr.T,b_arr.T, p_arr.T, q_arr.T))


@mem
def f_lkFiqn(l, k, i, n):
    '''
    This function computes the values of the easy basis f_lk(F_i (q_n)) 
        (level 1 of SG) based on the Splines on Fractals paper.

    Args:
        l, k, i, n: correspond to the indices mentioned in the preceding 
            paragraph.

    Returns:
        values of the easy basis f_lk(F_i (q_n))
    '''
    # Finds the values of pi and qi from the recursion
    arr = big_recursion(l)
    p1 = copy.deepcopy(arr[2, l])
    q1 = copy.deepcopy(arr[3, l])
    
    # Find the values of easy basis from (5.2) of splines paper
    if i == n and i == k:
        return int(l == 0)
    if i == k and i != n:
        return gm.mpq(p1,(5**l))
    if i == n and i != k:
        return 0
    if k == n and k != i:
        return gm.mpq(p1,(5**l))

    return gm.mpq(q1,(5**l))


@mem
def f_jk(addr, j, k):
    '''
    This function calculates the value of the easy basis f_jk at a point 
        on SG addressed by addr. This code is based on the Splines on 
        Fractals paper.

    Args:
        j, k: indices mentioned in preceding paragraph.
        addr: address of evaluation point F_w(q_i) given as a string of 
            digits whose first digit is i and whose following digits are 
            those in w.

    Returns:
        value of f_jk(F_w(q_i))
    '''
    # Base case 1: f_jk(q_i)
    if len(addr) == 1:
        if j != 0:
            return 0
        return int(int(addr[0]) == k)
    
    # Base case 2: f_jk(F_iq_n)
    if len(addr) == 2:
        n = int(addr[0])
        i = int(addr[1])
        return f_lkFiqn(j, k, i, n)
    
    # Inductive case, from (2.8) of splines paper
    last = int(addr[-1])
    addr = addr[:-1]
    outer_sum = 0
    for l in range(j+1):
        inner_sum = 0
        for n in range(3):
            inner_sum += f_lkFiqn(j-l, k, last, n)*f_jk(addr, l, n)
        inner_sum *= (gm.mpq(1, 5))**l
        outer_sum += inner_sum
    return outer_sum 


@mem
def p_jk(addr, j, k):
    '''
    This function calculates the value of the monomial basis p_jk at a 
        point on SG addressed by addr. This code is based on the Splines 
        on Fractals paper and Calculus on SG I.

    Args:
        j, k: indices mentioned in preceding paragraph.
        addr: address of evaluation point F_w(q_i) given as a string of 
            digits whose first digit is i and whose following digits are 
            those in w.

    Returns:
        value of p_jk(F_w(q_i))
    '''
    # Based on (2.14), (2.20) of the Calculus paper
    if k == 1:
        res = f_jk(addr, j, 0)
        for l in range(j+1):
            res += alpha(j-l)*(f_jk(addr, l, 1) + f_jk(addr, l, 2))
    if k == 2:
        res = 0
        for l in range(j+1):
            res += beta(j-l)*(f_jk(addr, l, 1) + f_jk(addr, l, 2))
    if k == 3:
        res = 0
        for l in range(j+1):
            res += gamma(j-l)*(f_jk(addr, l, 1) - f_jk(addr, l, 2))

    return res


def generate_W(level, deg, frac=True):
    '''
    This function calculates the values of the easy basis f_jk up to a certain 
        degree at a given level of SG.

    Args:
        level: level of SG required
        deg: maximum degree monomial required

    Returns:
        W: np.array with dimensions 3 x 3^(level+1) x deg+1.
        The kth page of W has values of the easy basis f_jk (j = 0...deg) 
            at each of the 3^(level+1) TLR indices of points on the given level 
            of SG.
    '''
    # Computes big_recursion at the highest degree, and lower values are 
    #   stored for later use
    print('Computing Big_Recursion... this may take some time')
    big_recursion(deg+1)

    # Initialize the W array based on how things should be stored.
    if frac:
        W = np.empty((3, 3**(level + 1), deg+1), dtype='object')
    else:
        W = np.zeros((3, 3**(level + 1), deg+1))
    
    # Main loop to fill in the values of the W array
    for i in tqdm.tqdm(range(3**(level + 1))):
        # This preparation is due to the different address structure 
        #   used in this file an in util.py
        addr = address_from_index(level, i+1)
        addr = np.flip(addr)
        addr = ''.join(str(int(x)) for x in addr)
        for j in tqdm.tqdm(range(deg + 1)):
            for k in tqdm.tqdm(range(1, 4)):
                W[k-1, i, j] = f_jk(addr, j, k-1)
    return W


def generate_T(level, deg, frac=True):
    '''
    This function calculates the values of the monomials up to a certain 
        degree at a given level of SG.

    Args:
        level: level of SG required
        deg: maximum degree monomial required
        frac: Boolean representing whether the coefficients should 
            remain as fractions or should be converted to floating point 
            numbers at the end of all calculations.

    Returns:
        T: np.array with dimensions 3 x 3^(level+1) x deg+1.
            The kth page of T has values of the monomials P_jk 
            (j = 0...deg) at each of the 3^(level+1) TLR indices of 
            points on the given level of SG.
    '''

    # Computes big_recursion at the highest degree, and lower values 
    #   are stored for later use.
    print('Computing Big_Recursion... this may take some time')
    big_recursion(deg+1)
    
    # Initialize the T array based on how things should be stored.
    if frac:
        T = np.empty((3, 3**(level + 1), deg+1), dtype='object')
    else:
        T = np.zeros((3, 3**(level + 1), deg+1))

    # Main loop to fill in the values of the T array
    for i in tqdm.tqdm(range(3**(level + 1))):
        # This preparation is due to the different address structure 
        #   used in this file an in util.py
        addr = address_from_index(level, i+1)
        addr = np.flip(addr)
        addr = ''.join(str(int(x)) for x in addr)
        for j in tqdm.tqdm(range(deg + 1)):
            for k in tqdm.tqdm(range(1, 4)):
                T[k-1, i, j] = p_jk(addr, j, k)
    return T


# def generate_T_loop(i, j, k, level):
#     addr = address_from_index(level, i+1)
#     addr = np.flip(addr)
#     addr = ''.join(str(int(x)) for x in addr)
#     return p_jk(addr, j, k)

# def generate_T_parallel(level, deg, frac=True, cores=4):
#     dtype = 'object' if frac else np.float64 
#     T = Parallel(n_jobs=cores)(delayed(generate_T_loop)(i=i, j=j, k=k, level=level) for i in range(3**(level + 1)) for j in range(deg +1) for k in range(1, 4))
#     return np.array(T, dtype=dtype).reshape(((3**(level + 1), deg+1, 3)))

# level = 5
# deg = 10
# frac = 0
# cores = 4
#start = time.time()

#T1 = generate_T(level, deg, frac)
# print('Non-parallel: ', time.time() - start)
# T2 = generate_T_parallel(level, deg, frac, cores=cores)
# print('Parallel (' + str(cores) + ' cores): ' , time.time() - start)
# #print(T2)
# #print(p_jk('10000000', 1, 1))
