import numpy as np
from recursions import alpha, beta, gamma, eta, ap
import gmpy2 as gm
import itertools
from monomials import p_jk, generate_T, f_jk
from Polynomial import Polynomial
import sympy as sp
from util import HiddenPrints, bmatrix
from matplotlib import pyplot as plt
import tqdm

'''
This file contains functions used to study higher order quadrature on SG using n-harmonic splines
'''

def make_vortex(max_level):
    '''
        This function generates a list of addresses for the vertices chosen for
        interpolation using the "vortex method."

        Args:
            max_level: integer representing the level of SG used for interpolation
        
        Returns:
            list of addresses for the vertices chosen for
            interpolation using the "vortex method.
    '''
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
    '''
        This function creates the interpolation matrix used for determining
        quadrature weights using up to j-harmonic splines.

        Args:
            j: maximum degree of harmonic spline used
        
        Returns:
            np.array that is the interpolation matrix used for determining
            quadrature weights using up to j-harmonic splines.
    '''
    addresses = make_vortex(j)

    q_mat = np.zeros((3*j+3, 3*j+3))
    
    for ind1 in range(3*j+3):
        for ind2 in range(3*j+3):
            jj = int(np.floor(ind1/3)) 
            kk = int(ind1 % 3 + 1)
            with HiddenPrints():
                q_mat[ind1, ind2] = p_jk(addresses[ind2], jj, kk)
    return q_mat

def ints_vec(j):
    '''
        This function computes the integrals of {P_01, P_02, P_03, ..., P_j2, P_j3}.

        Args:
            j: maximum degree of polynomial
        
        Returns:
            np.array of integrals of {P_01, P_02, P_03, ..., P_j2, P_j3}


    '''
    res = np.zeros(3*j+3)
    for ind in range(3*j+3):
        jj = int(np.floor(ind/3)) 
        kk = int(ind % 3 + 1)
        res[ind] = Polynomial.basis_inner(jj, kk, 0, 1)
    return res

def get_weights(j):
    '''
        This function determines the quadrature weights for n-Harmonic splines quadrature
        using the P_jk basis

        Args:
            j: maximum degree of polynomial used for quadrature

        Returns:
            np.array of quadrature weights for the basis {P_01, P_02, P_03, ..., P_j2, P_j3}

    '''
    return np.linalg.solve(quad_mat(j), ints_vec(j))

#print(np.linalg.inv(quad_mat(5)))
#print(make_vortex(3))
#print(sum(get_weights(20)))


def quad_int(func, j):
    '''
        This function calculates the quadrature integral for a given function.

        Args:
            func: function whose input is the address of a point on SG 
            (according to the convention in monomials.py)
            j: maximum degree of polynomial used for quadrature

        Returns:
            approximate integral of func

    '''
    addresses = make_vortex(j)
    func_vals = np.zeros(3*j+3)
    for i in range(len(addresses)):
        func_vals[i] = func(addresses[i])
    weights = get_weights(j)
    return func_vals.dot(weights)



def make_block(ind1, ind2):
    '''
        This function is used to generate the blocks in the block-matrix form 
        of the generalized vortex method for interpolation.

        Args:
            ind1, ind2: block matrix indices in the block form of the 
            generalized vortex method for interpolation

        Returns:
            sp.Matrix block C_ij in the block form of the 
            generalized vortex method for interpolation

    '''

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
    '''
        This function creates the block form of the generalized vortex method for interpolation

        Args:
            j: maximum degree of polynomial used for quadrature

        Returns:
            sp.Block_Matrix representing the generalized vortex method for interpolation

    '''
    big_list = []
    for ind1 in range(1,j+1):
        small_list = []
        for ind2 in range(1,j+1):
            small_list.append(make_block(ind1, ind2))
        big_list.append(small_list)
    
    return sp.BlockMatrix(big_list)

def block_to_regular(mat):
    '''
        This function converts an sp.Block_Matrix to an sp.Matrix
        with the same elements

        Args:
            mat: sp.Block_Matrix to be converted

        Returns:
            sp.Matrix with the same elements as mat
    '''
    res = sp.zeros(mat.rows, mat.cols)
    for i in range(mat.rows):
        for j in range(mat.cols):
            res[i, j] = mat[i, j]
    return res


