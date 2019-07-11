import numpy as np
from recursions import alpha, beta, gamma, eta, ap
import gmpy2 as gm
import itertools
from monomials import p_jk, generate_T
from Polynomial import Polynomial


def make_vortex(max_level):
    addr0 = []
    addr1 = []
    addr2 = []
    for i in range(max_level+1):
        ad0 = '0'
        ad1 = '1'
        ad2 = '2'
        for _ in itertools.repeat(None, i):
            ad0 += '2'
            ad1 += '0'
            ad2 += '1'
        addr0.append(ad0)
        addr1.append(ad1)
        addr2.append(ad2) 
    return addr0 + addr1 + addr2

def quad_mat(j):
    addresses = make_vortex(j)

    q_mat = np.zeros((3*j+3, 3*j+3))
    
    for ind1 in range(3*j+3):
        for ind2 in range(3*j+3):
            jj = int(np.floor(ind1/3)) 
            kk = int(ind1 % 3 + 1)
            q_mat[ind1, ind2] = p_jk(addresses[ind2], jj, kk)
    return q_mat

def ints_vec(j):
    res = np.zeros(3*j+3)
    for ind in range(3*j+3):
        jj = int(np.floor(ind/3)) 
        kk = int(ind % 3 + 1)
        res[ind] = Polynomial.basis_inner(jj, kk, 0, 1)
    return res
def get_weights(j):
    return np.linalg.solve(quad_mat(j), ints_vec(j))

#print(np.linalg.inv(quad_mat(5)))
#print(make_vortex(3))
#print(sum(get_weights(20)))


def quad_int(func, j):
    addresses = make_vortex(j)
    func_vals = np.zeros(3*j+3)
    for i in range(len(addresses)):
        func_vals[i] = func(addresses[i])
    weights = get_weights(j)
    return func_vals.dot(weights)

# p_52 = lambda add: p_jk(add, 5, 2)

# quad = quad_int(p_52, 4)
# actual = Polynomial.basis_inner(5, 2, 0, 1)
# print(quad)
# print(float(actual))
# print(abs(quad-actual))



