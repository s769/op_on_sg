import numpy as np
from recursions import alpha, beta, gamma, eta, ap
import gmpy2 as gm
import itertools
from monomials import p_jk, generate_T, f_jk
from Polynomial import Polynomial
import sympy as sp


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

def make_block(ind1, ind2):


    # ad0, ad1, ad2 = '0', '1', '2'
    ad0 = '0'
    for _ in itertools.repeat(None, ind1):
        ad0 += '1'
        # ad1 += '2'
        # ad2 += '0'

    
    a = f_jk(ad0, ind2, 1)
    b = f_jk(ad0, ind2, 2)
    c = f_jk(ad0, ind2, 3)
    return sp.Matrix([[a, b, c], [c, a, b], [b, c, a]])

def make_big_mat(j):
    big_list = []
    for ind1 in range(1,j+1):
        small_list = []
        for ind2 in range(1,j+1):
            small_list.append(make_block(ind1, ind2))
        big_list.append(small_list)
    
    return sp.BlockMatrix(big_list)
def block_to_regular(mat):
    res = sp.zeros(mat.rows, mat.cols)
    for i in range(mat.rows):
        for j in range(mat.cols):
            res[i, j] = mat[i, j]
    return res
res = block_to_regular(make_big_mat(2))

# for i in range(10):
#     print(sp.block_collapse(res**i))

# #print(res.eigenvects())
# print(res.eigenvals())

# C1 = make_block(3, 2)
# C2 = make_block(2, 3)
# print(C1*C2-C2*C1)
#print(sp.block_collapse(res.inverse()))

