import numpy as np
from Polynomial import Polynomial
from innerprods import lis2str
import scipy.io
import sympy as sp
from sympy import Rational as Rat
import tqdm

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

    Returns:
        np.array of coefficients of the orthogonal polynomials with 
            respect to the basis {P_0k, P_1k,..., P_nk}. Each row in 
            this array is a polynomial, and there are n+1 rows and n+1 
            columns.
            If normalized is True, the polynomials will be normalized. 
            Otherwise, the polynomials will be monic.
    '''
    print('Building Gram Matrix ... this may take some time')
    Polynomial.build_condensed_GM(n+1, k, lam)
#   basis_mat = np.zeros((n+1, 3*n+3))

#   for i in range(n+1):
#     basis_mat[i, 3*i + k - 1] = 1
    basis_mat = sp.eye(n+1)
    #o_basis_mat = np.zeros((n+1, 3*n+3))
    o_basis_mat = sp.zeros(n+1, n+1)
    

    o_basis_mat[0,:] = basis_mat[0,:]
    GM = Polynomial.GM[lis2str(lam)][:n+1, :n+1]
    print('Orthogonalizing Using Gram-Schmidt')
    for r in tqdm.tqdm(range(1, n+1)):
        #res = np.zeros(n+1)

        u_r = basis_mat[r,:]
        for i in range(r):
            v_i = o_basis_mat[i,:]
            proj = Polynomial.fast_inner(u_r.T, v_i.T, GM)
            norm = Polynomial.fast_inner(v_i.T, v_i.T, GM)
            u_r -= (proj/norm)*v_i
        o_basis_mat[r,:] = u_r#basis_mat[r] - res

    if normalized:
        for i in range(n+1):
            norm = Polynomial.fast_inner(o_basis_mat[i,:].T, o_basis_mat[i,:].T,
                                         GM)
            o_basis_mat[i,:] = o_basis_mat[i,:]/sp.sqrt(norm[0])
    return o_basis_mat  if frac else np.array(o_basis_mat).astype(np.float64)# , o_basis_mat[:, k-1::3]


# j = 5
# k = 3
# normalized = 1
# ops_sob = generate_op(j, k, normalized, lam=np.array([1]), frac=0)
# # ops_leg = generate_op(j, k, normalized, lam=np.array([0]))
# print(ops_sob)
# # #arr1 = ops[2]
# # # arr2 = ops[3]
# # # Polynomial.fast_inner(arr1, arr2, Polynomial.GM)]
# # #print(ops)


# # poly_arr = np.array([Polynomial(r, j, k) for r in ops[0]])

# # for k in range(j):
# #   sn1 = poly_arr[k+1].value()
# #   sn = poly_arr[k].value()
# #   dsn1 = poly_arr[k+1].dnvalue()
# #   dsn = poly_arr[k].dnvalue()
# #   print(sn1*dsn - sn*dsn1)
# #print(repr(ops))

# ops = np.zeros((2*j+2, j+1))
# for i in range(2*j+2):
#   if not i%2:
#     ops[i] = ops_leg[int(i/2)]
#   elif i%2:
#     ops[i] = ops_sob[int(i/2)]


#scipy.io.savemat('../data/coefs.mat', dict(ops=ops_leg))
# scipy.io.savemat('coefs2.mat', dict(ops=ops_leg))
# scipy.io.savemat('coefscomb.mat', dict(ops=ops))


# ops_sob

# ops = generate_op(j, k, 1, lam=np.ones(C))
# ops
# ops_sob
