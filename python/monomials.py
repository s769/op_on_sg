import numpy as np
from recursions import mem, mem2, alpha, beta, gamma, zeros, eye
from util import address_from_index, progress
import sympy as sp
from sympy import Rational as Rat
import tqdm
#from joblib import Parallel, delayed
import time
import gmpy2 as gm
import copy
#import const
'''
This file contains the functions that will compute the values of the 
    monomial basis on a given level of SG. The recursive functions are 
    memoized (see recursions.py).
'''


@mem2
def big_recursion(j):

    '''
    This function computes the coefficients a_j, b_j, p_j, q_j found in 
        the Splines on Fractals paper.

    Args:
        j: index of coefficients

    Returns:
        coefficients a_j, b_j, p_j, q_j
    '''
    p_arr = zeros(j+1,1)
    q_arr = zeros(j+1,1)
    a_arr = zeros(j+1,1)
    b_arr = zeros(j+1,1)

    a_arr[0] = gm.mpq(7, 45)
    b_arr[0] = gm.mpq(4, 45)
    p_arr[0] = gm.mpq(2, 5)
    q_arr[0] = gm.mpq(1, 5)

    if j == 0:
        return np.vstack((a_arr, b_arr, p_arr, q_arr))
        # return np.array([gm.mpq(7,45), gm.mpq(4,45), gm.mpq(2,5), gm.mpq(1,5)])

    

    
    # vec2 = zeros(2, 1)
    for l in range(1, j+1):
        res3 = 0
        res4 = 0
        for k in range(l):
            p = p_arr[k]
            q = q_arr[k]
            a = a_arr[l-k-1]
            b = b_arr[l-k-1]
            # p, q = big_recursion(l)[-2:]
            # a, b = big_recursion(j-1-l)[:2]
            res3 += (4*a + 3*b)*p + (a + 2*b)*q
            res4 += (2*a + 4*b)*p + (3*a + b)*q
        
    # b = big_recursion(j-1)[1]
        b = b_arr[l-1]
        p_arr[l] = -gm.mpq(2,5)*b - gm.mpq(1,5)*res3
        q_arr[l] = -gm.mpq(1,5)*b - gm.mpq(1,5)*res4

        res1 = 0
        res2 = 0
        vec = zeros(2,1)
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
    

    # for l in range(j):
    #     res1 = 0
    #     res2 = 0
    #     vec = zeros(2,1)
    #     if l == 0:
    #         p, q = p_arr[-1], q_arr[-1]
    #     else:
    #         p = p_arr[j-l]
    #         q = q_arr[j-l]
            
    #         # p, q = big_recursion(j-l)[-2:]
    #     a = a_arr[l]
    #     b = b_arr[l]
    #     # a, b = big_recursion(l)[:2]
    #     res1 += (2*p + q)*(a+2*b)
    #     res2 += (p + 2*q)*(a + 2*b)
    # vec[0] = gm.mpq(2,(3*(5**j - 1)))*res1
    # vec[1] = gm.mpq(10,(3*(5**(j+1) - 1)))*res2
    # coefs_inv = np.array([[gm.mpq(7,15), gm.mpq(-2,15)], [gm.mpq(4,15), gm.mpq(1,15)]])
    # ab_arr = coefs_inv@vec
    # # ab_arr = np.append(ab_arr, vec2)
    # a_arr[-1] = ab_arr[0]
    # b_arr[-1] = ab_arr[1]
  
    
    return np.vstack((a_arr.T,b_arr.T, p_arr.T, q_arr.T))

# rec_arr = big_recursion_master(j_max+1)
# def big_recursion(j):
#     return rec_arr#[:,:j+1] if j < j_max else big_recursion_master(j)
    

#arr = big_recursion(3)
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

    # print(arr)
    # print(type(arr))
    arr = big_recursion(l)
    p1 = copy.deepcopy(arr[2, l])
    q1 = copy.deepcopy(arr[3, l])
    
    # p, q = big_recursion(l)[-2:]
    #print(p1)
    #print(q1)
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
    #print('before', addr)
    if len(addr) == 1:
        if j != 0:
            return 0
        return int(int(addr[0]) == k)
    if len(addr) == 2:
        n = int(addr[0])
        i = int(addr[1])
        return f_lkFiqn(j, k, i, n)
    last = int(addr[-1])
    addr = addr[:-1]
    resouter = 0
    #print('after', addr) 
    for l in range(j+1):
        resinner = 0
        for n in range(3):
            # print('first' , f_lkFiqn(j-l, k, last, n)*f_jk(addr, l, n))
            # print('resinner after first', resinner)
            # print(f_lkFiqn(j-l, k, last, n)*f_jk(addr, l, n))
            # print(resinner + f_lkFiqn(j-l, k, last, n)*f_jk(addr, l, n))
            #ad = addr
            a = f_lkFiqn(j-l, k, last, n)*f_jk(addr, l, n)
            #print('a', ad, l, n, j-l, k, last,  a)
            resinner += a#f_lkFiqn(j-l, k, last, n)*f_jk(addr, l, n)

            #print('inner res', resinner)
        resinner *= (gm.mpq(1, 5))**l
        #print('outer res' , resinner)
        resouter += resinner
        #print('outer outer res', resouter)
    return resouter 
#print(f_jk('01221', 2, 3))
#print(big_recursion(100))
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

# ans = f_lkFiqn(0, 1, 1, 0)*f_lkFiqn(0,0,1,0)+f_lkFiqn(0,1,1,1)*f_lkFiqn(0,1,1,0)+f_lkFiqn(0,1,1,2)*f_lkFiqn(0,2,1,0)
# print(ans)
#print(f_lkFiqn(0,2,1,0))
#print(f_jk('0', 0, 3))


def generate_T(level, deg, frac=True):
    '''
    This function calculates the values of the monomials up to a certain 
        degree at a given level of SG.

    Args:
        level: level of SG required
        deg: maximum degree monomial required

    Returns:
        T: np.array with dimensions 3 x 3^(level+1) x deg+1.
        The kth page of T has values of the monomials P_jk (j = 0...deg) 
            at each of the 3^(level+1) TLR indices of points on the given level 
            of SG.
    '''
    #computes big_recursion at the highest degree so lower values are stored for later use
    # print('Computing Big_Recursion... this may take some time')
    big_recursion(deg+1)
    if frac:
        #T = np.empty((3**(level + 1), deg+1, 3), dtype='object')
        T = np.empty((3, 3**(level + 1), deg+1), dtype='object')
    else:
        #T = np.zeros((3**(level + 1), deg+1, 3))
        T = np.zeros((3, 3**(level + 1), deg+1))
    for i in tqdm.tqdm(range(3**(level + 1))):
        # This preparation is due to the different address structure 
        #   used in this file an in util.py
        addr = address_from_index(level, i+1)
        addr = np.flip(addr)
        addr = ''.join(str(int(x)) for x in addr)
        for j in tqdm.tqdm(range(deg + 1)):
            for k in tqdm.tqdm(range(1, 4)):
                #T[i, j, k-1] = p_jk(addr, j, k)
                T[k-1,i, j] = p_jk(addr, j, k)
        #progress(i, 3**(level+1), status='computing monomial values')
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
