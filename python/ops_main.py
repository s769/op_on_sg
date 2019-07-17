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





def leg_ops_recursion(j, k, frac=True):
    '''
        This function uses the three term recursion from the Kasso Tuley paper to generate the first j
        Legendre orthogonal polynomials.

        Args:
            j: maximum degree of polynomials
            k: family of monomials to use in the construction of the orthogonal polynomials 
            (only k = 2,3 supported currently)
            frac: Boolean representing whether the coefficients should remain as fractions or should be
            converted to floating point numbers at the end of all calculations

        Returns:
            np.array of coefficients of the Legendre orthogonal polynomials with 
                respect to the basis {P_0k, P_1k,..., P_jk}. Each row in 
                this array is a polynomial, and there are j+1 rows and j+1 
                columns.
    '''

    if k == 1:
        print('This method is currently only available for k = 2 or 3.')
        return
    # this is so the indices match    
    j -= 1
    dtype = object if frac else np.float64
    o_basis_mat = np.empty((j+2, j+2), dtype=dtype)

    first_mat = generate_op(1, k, normalized=False, lam=np.array([0]), frac=frac)
    const = gm.mpz(0) if frac else 0
    first_mat = np.pad(first_mat, ((0,0), (0, j)),'constant', constant_values=(const,))
    o_basis_mat[:2] = first_mat

    func_array = gamma_array if k == 3 else beta_array

    func_arr = func_array(j+3)
    Polynomial.build_condensed_GM(j+2, k, np.array([0]))
    GM = Polynomial.GM[lis2str(np.array([0]))][:j+2, :j+2]

    for ind in range(1,j+1):
        func_vec = func_arr[1:ind+2]
        omega_vec = o_basis_mat[ind, :ind+1]
        zeta_ind = -1/func_arr[0]*func_vec.dot(omega_vec)
        f_ind = np.insert(omega_vec, 0, zeta_ind)
        d_ind2 = gm.mpq(1, Polynomial.fast_inner(o_basis_mat[ind,:ind+1], o_basis_mat[ind,:ind+1], GM[:ind+1, :ind+1]))
        d_indm2 = gm.mpq(1, Polynomial.fast_inner(o_basis_mat[ind-1,:ind], o_basis_mat[ind-1,:ind], GM[:ind, :ind]))
        b_ind = d_ind2*Polynomial.fast_inner(f_ind, o_basis_mat[ind,:ind+2], GM[:ind+2, :ind+2])
        c_ind = gm.mpq(d_indm2, d_ind2)

        new_vec = f_ind - b_ind*o_basis_mat[ind, :ind+2] - c_ind*o_basis_mat[ind-1, :ind+2]

        o_basis_mat[ind+1] = np.pad(new_vec, (0, j-ind), 'constant', constant_values=(const,))

    return o_basis_mat


def sob_ops_recursion(j, k, frac=True, leg_omegas=None):
    '''
        This function uses the three term recursion we developed to generate the first j
        Sobolev orthogonal polynomials.

        Args:
            j: maximum degree of polynomials
            k: family of monomials to use in the construction of the orthogonal polynomials 
            (only k = 2,3 supported currently)
            frac: Boolean representing whether the coefficients should remain as fractions or should be
            converted to floating point numbers at the end of all calculations

        Returns:
            np.array of coefficients of the Sobolev orthogonal polynomials with 
                respect to the basis {P_0k, P_1k,..., P_jk}. Each row in 
                this array is a polynomial, and there are j+1 rows and j+1 
                columns.
    '''
    if k == 1:
        print('This method is currently only available for k = 2 or 3.')
        return
    # this is so the indices match    
    j -= 1
    dtype = object if frac else np.float64
    o_basis_mat = np.empty((j+2, j+2), dtype=dtype)

    first_mat = generate_op(1, k, normalized=False, frac=frac)
    const = gm.mpz(0) if frac else 0
    first_mat = np.pad(first_mat, ((0,0), (0, j)),'constant', constant_values=(const,))
    o_basis_mat[:2] = first_mat

    
    func_array = gamma_array if k == 3 else beta_array

    func_arr = func_array(j+3)
    Polynomial.build_condensed_GM(j+2, k, np.array([1]))
    GM = Polynomial.GM[lis2str(np.array([1]))][:j+2, :j+2]

    if leg_omegas is None:
        leg_omegas = leg_ops_recursion(j+1, k, frac=frac)
    elif isinstance(leg_omegas, tuple):
        filename, arr = leg_omegas
        leg_omegas = np.load(filename, allow_pickle=frac)[arr]
    
    for ind in range(1,j+1):
        func_vec = func_arr[1:ind+2]
        omega_vec = leg_omegas[ind, :ind+1]
        zeta_ind = -1/func_arr[0]*func_vec.dot(omega_vec)
        f_ind = np.insert(omega_vec, 0, zeta_ind)
        a_ind = Polynomial.fast_inner(f_ind, o_basis_mat[ind,:ind+2], GM[:ind+2, :ind+2])
        b_ind = Polynomial.fast_inner(f_ind, o_basis_mat[ind-1,:ind+2], GM[:ind+2, :ind+2])
        a_ind = gm.mpq(a_ind, Polynomial.fast_inner(o_basis_mat[ind, :ind+1], o_basis_mat[ind,:ind+1], GM[:ind+1, :ind+1]))
        b_ind = gm.mpq(b_ind, Polynomial.fast_inner(o_basis_mat[ind-1, :ind], o_basis_mat[ind-1,:ind], GM[:ind, :ind]))
        new_vec = f_ind - a_ind*o_basis_mat[ind, :ind+2] - b_ind*o_basis_mat[ind-1, :ind+2]

        o_basis_mat[ind+1] = np.pad(new_vec, (0, j-ind), 'constant', constant_values=(const,))

    return o_basis_mat

