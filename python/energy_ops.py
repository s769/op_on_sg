import numpy as np
import gmpy2 as gm
import tqdm

from recursions import alpha_array, beta_array, gamma_array, eta_array,\
    ap_array, zeros_gm, eye_gm, alpha, beta, gamma, eta
from innerprods import lis2str, symmetrize
from Polynomial import Polynomial
from util import HiddenPrints
from monomials import p_jk, norm_p_jk
import sys, os


# def p_jkqi(j, k, i):
#     if i == 0:
#         if k == 2 or k == 3: return gm.mpz(0)
#         return j == 0
#     if k == 1:
#         return alpha(j)
#     if k == 2:
#         return beta(j)
#     if k == 3:
#         return gamma(j) if i == 1 else -gamma(j)

# def dn_p_jkqi(j, k, i):
#     if i == 0:
#         if k == 1 or k == 3: return gm.mpz(0)
#         return j == 0
#     if k == 1:
#         return eta(j)
#     if k == 2:
#         return -alpha(j)
#     if k == 3:
#         return 3*eta(j+1) if i == 1 else -3*eta(j+1)

def energy_inner_pjkpmn(j, k, m, n, lam=np.array([1])):
    '''
        Calculates energy inner product between monomials P_{j,k} and P_{m, n}.
        Uses the Gauss-Green formula so that 
            <P_{j,k}, P_{m,n}>_energy = <P_{j,k}, P_{m,n}>_L2 - lam*<P_{j-1,k}, P_{m,n}>_L2
                                                    + lam*sum_V0{P_{m,n}*d_n(P_{j,k})}
        Args:
            j, k, m, n: Represent monomial indices P_{j,k} and P_{m, n}
            lam: represents weight given to the energy part of the inner product
        Returns: 
            <P_{j,k}, P_{m,n}>_energy
    '''
    lam = lam[0]
    res = Polynomial.basis_inner(j, k, m, n, lam=np.array([0]))
    res -= lam*Polynomial.basis_inner(j-1, k, m, n, lam=np.array([0]))
    for i in range(2):
        with HiddenPrints():
            res += lam*norm_p_jk(str(i), j, k)*p_jk(str(i), m, n)
    return res

def generate_energy_ops(n, k, normalized=False, lam = np.array([1]), frac=True):
    '''
    Generates orthogonal polynomials with respect to the energy inner product. The Gram-Schmidt algorithm 
        is implemented here.

    Args:
        n: Maximum degree of orthogonal polynomial.
        k: family of monomials to use in Gram-Schmidt (k = 1, 2, or 3)
        normalized: Boolean representing whether the resulting polynomials 
            should be normalized or monic.
        lam: np.array of lambda values for the generalized Sobolev inner 
            product. The default value is 1 (corresponding to the regular 
            Sobolev inner product). If lam = np.array([0]), 
            this is the L2 inner product.
        frac: Boolean representing whether the coefficients should remain as fractions or should be
        converted to floating point numbers at the end of all calculations.

    Returns:
        np.array of coefficients of the orthogonal polynomials with 
            respect to the basis {P_0k, P_1k,..., P_nk}. Each row in 
            this array is a polynomial, and there are n+1 rows and n+1 
            columns.
            If normalized is True, the polynomials will be normalized. 
            Otherwise, the polynomials will be monic. If normalized is True, frac must be False
            to obtain normalized coefficients.
    
    '''

    print('Building Gram Matrix ... this may take some time')
    EGM = zeros_gm(n+1, n+1)

    for ind1 in range(n+1):
        for ind2 in range(n+1):
            if ind1 <= ind2:
                EGM[ind1, ind2] = energy_inner_pjkpmn(ind1, k, ind2, k, lam=lam)

    EGM = symmetrize(EGM)


    basis_mat = eye_gm(n+1)
    o_basis_mat = zeros_gm(n+1, n+1)

    o_basis_mat[0] = basis_mat[0]

    print('Orthogonalizing Using Gram-Schmidt')
    for r in tqdm.tqdm(range(1, n+1)):
        u_r = basis_mat[r]
        for i in range(r):
            v_i = o_basis_mat[i]
        
            proj = Polynomial.fast_inner(u_r, v_i, EGM)
            norm = Polynomial.fast_inner(v_i, v_i, EGM)
            u_r -= (proj/norm)*v_i
        o_basis_mat[r] = u_r

    if frac and normalized:
        print('Normalization requires conversion to float. Please set frac = False.')
        print('Generating non-normalized coefficients now...')
    if not frac:
        if normalized:
            o_basis_arr = np.zeros((n+1, n+1))
            print('Normalizing')
            for i in tqdm.tqdm(range(n+1)):
                norm = Polynomial.fast_inner(o_basis_mat[i], o_basis_mat[i],
                                        EGM)
                o_basis_arr[i] = o_basis_mat[i]/gm.sqrt(norm)
            return o_basis_arr
        return np.array(o_basis_mat, dtype=np.float64)

    return o_basis_mat

