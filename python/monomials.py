import numpy as np
import progressbar
from recursions import mem, alpha, beta, gamma
from util import address_from_index

'''
This file contains the functions that will compute the values of the monomial basis on a given level of SG.
The recursive functions are memoized (see recursions.py).
'''


@mem
def big_recursion(j):
    '''
    This function computes the coefficients a_j, b_j, p_j, q_j found in the Splines on Fractals paper.

    Args:
        j: index of coefficients

    Returns:
        coefficients a_j, b_j, p_j, q_j
    '''
    if j == 0:
        return np.array([7/45, 4/45, 2/5, 1/5])

    res3 = 0
    res4 = 0
    vec2 = np.zeros(2)
    for l in range(j):
        p, q = big_recursion(l)[-2:]
        a, b = big_recursion(j-1-l)[:2]
        res3 += (4*a + 3*b)*p + (a + 2*b)*q
        res4 += (2*a + 4*b)*p + (3*a + b)*q
    b = big_recursion(j-1)[1]
    vec2[0] = -2/5*b - 1/5*res3
    vec2[1] = -1/5*b - 1/5*res4

    vec = np.zeros(2)
    res1 = 0
    res2 = 0

    for l in range(j):
        if l == 0:
            p, q = vec2
        else:
            p, q = big_recursion(j-l)[-2:]
        a, b = big_recursion(l)[:2]
        res1 += (2*p + q)*(a+2*b)
        res2 += (p + 2*q)*(a + 2*b)
    vec[0] = 2/(3*(5**j - 1))*res1
    vec[1] = 10/(3*(5**(j+1) - 1))*res2
    coefs = np.array([[1, 2], [-4, 7]])
    ab_arr = np.linalg.inv(coefs)@vec

    return np.append(ab_arr, vec2)


def f_lkFiqn(l, k, i, n):
    '''
    This function computes the values of the easy basis f_lk(F_i (q_n)) (level 1 of SG) based on the Splines on Fractals
    paper.

    Args:
        l, k, i, n: correspond to the indices mentioned in the preceding paragraph.

    Returns:
        values of the easy basis f_lk(F_i (q_n))

    '''
    p, q = big_recursion(l)[-2:]

    if i == n and i == k:
        return int(l == 0)
    if i == k and i != n:
        return p/5**l
    if i == n and i != k:
        return 0
    if k == n and k != i:
        return p/5**l

    return q/5**l


@mem
def f_jk(addr, j, k):
    '''
    This function calculates the value of the easy basis f_jk at a point on SG addressed by addr. This code is based on
    the Splines on Fractals paper.

    Args:
        j, k: indices mentioned in preceding paragraph.
        addr: address of evaluation point F_w(q_i) given as a string of digits whose first digit is i and
        whose following digits are those in w.

    Returns:
        value of f_jk(F_w(q_i))

    '''

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
    res = 0

    for l in range(j+1):
        for n in range(3):
            res += f_lkFiqn(j-l, k, last, n)*f_jk(addr, l, n)
        res *= (1/5)**l

    return res


def p_jk(addr, j, k):
    '''
    This function calculates the value of the monomial basis p_jk at a point on SG addressed by addr. This code is based on
    the Splines on Fractals paper and Calculus on SG I.

    Args:
        j, k: indices mentioned in preceding paragraph.
        addr: address of evaluation point F_w(q_i) given as a string of digits whose first digit is i and
        whose following digits are those in w.

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
# #print(f_lkFiqn(0,2,1,0))
#print(f_jk('0', 0, 3))


def generate_T(level, deg):
    '''
    This function calculates the values of the monomials up to a certain degree at a given level of SG.

    Args:
        level: level of SG required
        deg: maximum degree monomial required

    Returns:
        T: np.array with dimensions 3 x 3^(level+1) x deg+1.
        The kth page of T has values of the monomials P_jk (j = 0...deg) at each of the 3^(level+1) TLR indices of points
        on the given level of SG.

    '''
    T = np.zeros((3, 3**(level + 1), deg+1))
    pbar = ProgressBar(widgets=[Percentage(), Bar()],
                       maxval=3**(level+1)+1).start()
    for i in range(3**(level + 1)):
        # This preparation is due to the different address structure used in this file an in util.py
        addr = address_from_index(level, i+1)
        addr = np.flip(addr)
        addr = ''.join(str(int(x)) for x in addr)
        for j in range(deg + 1):
            for k in range(1, 4):
                T[k-1, i, j] = p_jk(addr, j, k)
        pbar.update(i+1)
    pbar.finish()
    return T


#print(generate_T(7, 10))
print(p_jk('2', 0, 3))
