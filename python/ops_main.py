import numpy as np
import gmpy2 as gm
import tqdm

from recursions import alpha_array, beta_array, gamma_array, eta_array,\
    ap_array, zeros_gm, eye_gm
from innerprods import lis2str
from Polynomial import Polynomial


'''
This is the main file used for computing the orthogonal polynomials.
'''


def generate_op(n, k, normalized=True, lam=np.array([1]), frac=True):
    '''
    Generates orthogonal polynomials with respect to a generalized Sobolev 
        inner product. The Gram-Schmidt algorithm is implemented here.

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
    Polynomial.build_condensed_GM(n+1, k, lam)

    basis_mat = eye_gm(n+1)
    o_basis_mat = zeros_gm(n+1, n+1)
    

    o_basis_mat[0] = basis_mat[0]
    GM = Polynomial.GM[lis2str(lam)][:n+1, :n+1]
    print('Orthogonalizing Using Gram-Schmidt')
    for r in tqdm.tqdm(range(1, n+1)):

        u_r = basis_mat[r]
        for i in range(r):
            v_i = o_basis_mat[i]
            
            proj = Polynomial.fast_inner(u_r, v_i, GM)
            norm = Polynomial.fast_inner(v_i, v_i, GM)
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
                                            GM)
                o_basis_arr[i] = o_basis_mat[i]/gm.sqrt(norm)
            return o_basis_arr
        return np.array(o_basis_mat, dtype=np.float64)

    return o_basis_mat


