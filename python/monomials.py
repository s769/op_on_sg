import time
import copy

import numpy as np
import gmpy2 as gm
import tqdm

from recursions import mem, mem2, alpha, beta, gamma, eta, tau, zeros_gm
from util import address_from_index, rotate_address, index_from_address
import sys, os
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
    for l in range(1, j+1):
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



def f_jk_addr_list(addr, j, k):
    """
    Evaluate f_jk when the address provided is a list in the order
    [w1, w2, ..., wm, qn] for F_{w1}...F_{wm}q_{n}.

    Args:
        j, k: indices mentioned in preceding paragraph.
        addr: address of evaluation point F_w(q_i) given as a list of 
            indices of F and q as defined above.

    Returns:
        value of f_jk(F_w(q_i))
    """
    addr = np.flip(addr)
    addr = ''.join(str(int(x)) for x in addr)
    return(f_jk(addr, j, k))


def p_jk_addr_list(addr, j, k):
    """
    Evaluate p_jk when the address provided is a list in the order
    [w1, w2, ..., wm, qn] for F_{w1}...F_{wm}q_{n}.

    Args:
        j, k: indices mentioned in preceding paragraph.
        addr: address of evaluation point F_w(q_i) given as a list of 
            indices of F and q as defined above.

    Returns:
        value of p_jk(F_w(q_i))
    """
    addr = np.flip(addr)
    addr = ''.join(str(int(x)) for x in addr)
    return(p_jk(addr, j, k))


@mem
def norm_f_jk(addr, j, k):
    '''
    This function calculates the value of the normal derivative of the 
        easy basis partial_{n}f_jk at a point on SG addressed by addr. 

    Args:
        j, k: indices mentioned in preceding paragraph.
        addr: address of evaluation point F_w(q_n) given as a string of 
            digits whose first digit is i and whose following digits are 
            those in w.

    Returns:
        value of partial_{n}f_jk(F_w(q_n))
    '''

    # Finds the values of pi and qi from the recursion
    arr = big_recursion(j)
    a1 = copy.deepcopy(arr[0, j-1])
    b1 = copy.deepcopy(arr[1, j-1])

    # Base case: partial_{n}f_jk(q_i)
    if len(addr) == 1:
        if (j == 0) and (int(addr[0]) == k):
            return 2
        elif (j == 0) and (int(addr[0]) != k):
            return -1
        elif (j != 0) and (int(addr[0]) == k):
            return a1
        else:
            return b1
    
    # Inductive case
    last = int(addr[-1])
    addr = addr[:-1]
    outer_sum = 0
    for l in range(j+1):
        inner_sum = 0
        for n in range(3):
            inner_sum += f_lkFiqn(j-l, k, last, n)*norm_f_jk(addr, l, n)
        inner_sum *= (gm.mpq(1, 5))**l
        outer_sum += inner_sum
    outer_sum *= gm.mpq(5, 3)
    return outer_sum 


@mem
def norm_p_jk(addr, j, k):
    '''
    This function calculates the value of the normal derivative of the 
        monomial basis partial_{n}p_jk at a point on SG addressed by 
        addr.

    Args:
        j, k: indices mentioned in preceding paragraph.
        addr: address of evaluation point F_w(q_n) given as a string of 
            digits whose first digit is i and whose following digits are 
            those in w.

    Returns:
        value of partial_{n}p_jk(F_w(q_n))
    '''
    # Based on (2.14), (2.20) of the Calculus paper

    if k == 1:
        res = norm_f_jk(addr, j, 0)
        for l in range(j+1):
            res += alpha(j-l)*(norm_f_jk(addr, l, 1) + norm_f_jk(addr, l, 2))
    if k == 2:
        res = 0
        for l in range(j+1):
            res += beta(j-l)*(norm_f_jk(addr, l, 1) + norm_f_jk(addr, l, 2))
    if k == 3:
        res = 0
        for l in range(j+1):
            res += gamma(j-l)*(norm_f_jk(addr, l, 1) - norm_f_jk(addr, l, 2))

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
    for i in tqdm.tqdm(range(3**(level + 1)), file=sys.stdout):
        # This preparation is due to the different address structure 
        #   used in this file and in util.py
        addr = address_from_index(level, i+1)
        addr = np.flip(addr)
        addr = ''.join(str(int(x)) for x in addr)
        for j in tqdm.tqdm(range(deg + 1), file=sys.stdout):
            for k in tqdm.tqdm(range(1, 4), file=sys.stdout):
                W[k-1, i, j] = f_jk(addr, j, k-1)
    return W


def generate_norm_W(level, deg, frac=True):
    '''
    This function calculates the values of the normal derivative of easy 
        basis partial_{n}f_jk up to a certain degree at a given level of SG.

    Args:
        level: level of SG required
        deg: maximum degree monomial required

    Returns:
        norm_W: np.array with dimensions 3 x 3^(level+1) x deg+1.
            The kth page of W has values of the normal derivative of the 
            easy basis f_jk (j = 0...deg) at each of the 3^(level+1) TLR 
            indices of points on the given level of SG.
    '''
    # Computes big_recursion at the highest degree, and lower values are 
    #   stored for later use
    print('Computing Big_Recursion... this may take some time')
    big_recursion(deg+1)

    # Initialize the W array based on how things should be stored.
    if frac:
        norm_W = np.empty((3, 3**(level + 1), deg+1), dtype='object')
    else:
        norm_W = np.zeros((3, 3**(level + 1), deg+1))
    
    # Main loop to fill in the values of the norm_W array
    for i in tqdm.tqdm(range(3**(level + 1)), file=sys.stdout):
        # This preparation is due to the different address structure 
        #   used in this file an in util.py
        addr = address_from_index(level, i+1)
        addr = np.flip(addr)
        addr = ''.join(str(int(x)) for x in addr)
        for j in tqdm.tqdm(range(deg + 1), file=sys.stdout):
            for k in tqdm.tqdm(range(1, 4), file=sys.stdout):
                norm_W[k-1, i, j] = norm_f_jk(addr, j, k-1)
    return norm_W


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
    for i in tqdm.tqdm(range(3**(level + 1)), file=sys.stdout):
        # This preparation is due to the different address structure 
        #   used in this file and in util.py
        addr = address_from_index(level, i+1)
        addr = np.flip(addr)
        addr = ''.join(str(int(x)) for x in addr)
        for j in tqdm.tqdm(range(deg + 1), file=sys.stdout):
            for k in tqdm.tqdm(range(1, 4), file=sys.stdout):
                T[k-1, i, j] = p_jk(addr, j, k)
    return T


def generate_norm_T(level, deg, frac=True):
    '''
    This function calculates the values of the normal derivative of the 
        monomials up to a certain degree at a given level of SG.

    Args:
        level: level of SG required
        deg: maximum degree monomial required
        frac: Boolean representing whether the coefficients should 
            remain as fractions or should be converted to floating point 
            numbers at the end of all calculations.

    Returns:
        T: np.array with dimensions 3 x 3^(level+1) x deg+1.
            The kth page of T has values of the normal derivative of the 
            monomials partial_{n}P_jk (j = 0...deg) at each of the 
            3^(level+1) TLR indices of points on the given level of SG.
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
    for i in tqdm.tqdm(range(3**(level + 1)), file=sys.stdout):
        # This preparation is due to the different address structure 
        #   used in this file an in util.py
        addr = address_from_index(level, i+1)
        addr = np.flip(addr)
        addr = ''.join(str(int(x)) for x in addr)
        for j in tqdm.tqdm(range(deg + 1), file=sys.stdout):
            for k in tqdm.tqdm(range(1, 4), file=sys.stdout):
                T[k-1, i, j] = norm_p_jk(addr, j, k)
    return T


def generate_T_symmetric(level, deg, frac=True, T=None):
    '''
        This function calculates the values of the fully symmetric monomials up to a certain 
        degree at a given level of SG.

    Args:
        level: level of SG required
        deg: maximum degree monomial required
        frac: Boolean representing whether the coefficients should 
            remain as fractions or should be converted to floating point 
            numbers at the end of all calculations.
        T: array of monomial values at the required level and degree or 
        tuple of (filename of .npz/.npy file containing this array , array name key string)

    Returns:
        T: np.array with dimensions 3^(level+1) x deg+1 containing
            values of the monomials R_j = P_j1^(0) + P_j1^(1) + P_j1^(2)
            (j = 0...deg) at each of the 3^(level+1) TLR indices of 
            points on the given level of SG.
    '''

    if T is None:
        T = generate_T(level, deg, frac=frac)
    elif isinstance(T, (tuple)):
        filename, arr = T
        T = np.load(filename, allow_pickle=frac)[arr]
    
    dtype = object if frac else np.float64

    ST = np.empty(T.shape[1:], dtype=dtype)

    for row in range(ST.shape[0]):
        addr = address_from_index(level, row+1)
        addr1 = rotate_address(level, addr, 1)
        addr2 = rotate_address(level, addr, 2)
        row1 = index_from_address(level, addr1)
        row2 = index_from_address(level, addr2)
        ST[row] = T[0, row] + T[0, row1-1] + T[0, row2-1]


    return ST
