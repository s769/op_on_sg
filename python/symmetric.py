import numpy as np
import gmpy2 as gm
import tqdm

from recursions import alpha, beta, gamma, eta, zeros_gm, eye_gm
from innerprods import inner0_j1k1, symmetrize, lis2str
from Polynomial import Polynomial


# This file contains functions that generate coefficients for the 
# 	symmetric orthogonal polynomials with respect to  various 
# 	inner products.

def pnj1(nn, j, i):
  	'''
	This function calculates the value of P_{j,1}^(n) at q_i
		
	Args:
		nn: corresponds to the value of (n) in P_{j,1}^(n)
		j: corresponds to the value of j in P_{j,1}^(n)
		i: corresponds to the boundary index q_i
	
	Returns:
		Value of P_{j,1}^(n) at q_i
  	'''
  	return int(j == 0) if i == nn else alpha(j)

def d_pnj1(nn, j, i):
  	'''
	This function calculates the value of the normal derivative of 
		P_{j,1}^(n) at q_i.
	
	Args:
		nn: corresponds to the value of (n) in P_{j,1}^(n)
		j: corresponds to the value of j in P_{j,1}^(n)
		i: corresponds to the boundary index q_i
    
    Returns:
    	Value of the normal derivative of P_{j,1}^(n) at q_i
	'''
	return 0 if i == nn else eta(j)


def inner_pnj1pmk1(nn, j, m, k):
	'''
	This function calculates the value of the L2 inner product 
		<P_{j,1}^(n), P_{m, 1}^(k)>.
	
	Args:
		nn: corresponds to the value of (n) in P_{j,1}^(n)
		j: corresponds to the value of j in P_{j,1}^(n)
		m: corresponds to the value of m in P_{m, 1}^(k)
		k: corresponds to the value of k in P_{m, 1}^(k)
		
	Returns:
		Value of <P_{j,1}^(n), P_{m, 1}^(k)>
	'''
  	if nn == m: return inner0_j1k1(j,k)
  	res = 0
  	for l in range(j+1):
    	s1 = 0
    	for i in range(3):
      		s1 += pnj1(nn, j-l, i)*d_pnj1(m, k+l+1, i)\
        	-pnj1(m, k+l+1, i)*d_pnj1(nn, j-l, i)
    	res += s1
  	return res


def inner_rjrk(j, k, lam = np.array([1])):
  	'''
    Computes any Sobolev inner product of <R_j, R_k> (the basis functions 
		for the symmetric polynomials). 
		R_j = P_{j, 1}^(0) + P_{j, 1}^(1) + P_{j, 1}^(2)
        
	Args:
		j, k: correspond to j and k in R_j, R_k
			respectively.
		lam: np.array of lambda values for the generalized Sobolev 
			inner product. The default value is 1 (corresponding to 
			the regular Sobolev inner product). 
			If lam = np.array([0]), this is the L2 inner product.

	Returns:
		Sobolev inner product of <R_j, R_k> with given lambda 
			values.
  	'''
  	res = 0
  	for nn in range(3):
    	for m in range(3):
      		res = inner_pnj1pmk1(nn, j, m, k)
      		for i in range(1, len(lam) + 1):
        		res += lam[i-1]*inner_pnj1pmk1(nn, j-1, m, k-1)
  	return res


def generate_symm_ops(n, normalized=False, lam = np.array([1]), frac=True):
  	'''
    Generates symmetric orthogonal polynomials with respect to a 
		generalized Sobolev inner product. The Gram-Schmidt algorithm 
		is implemented here.

    Args:
        n: Maximum degree of orthogonal polynomial.
        normalized: Boolean representing whether the resulting polynomials 
            should be normalized or monic.
        lam: np.array of lambda values for the generalized Sobolev inner 
            product. The default value is 1 (corresponding to the regular 
            Sobolev inner product). If lam = np.array([0]), 
            this is the L2 inner product.
        frac: Boolean representing whether the coefficients should remain 
			as fractions or should be converted to floating point numbers 
			at the end of all calculations.

    Returns:
        np.array of coefficients of the orthogonal polynomials with 
            respect to the basis {R_0, R_1, ..., R_j}. Each row in 
            this array is a polynomial, and there are n+1 rows and n+1 
            columns.
            If normalized is True, the polynomials will be normalized.
            Otherwise, the polynomials will be monic. If normalized 
			is True, frac must be False to obtain normalized coefficients.
    
    '''
  	print('Building Gram Matrix ... this may take some time')
  	SGM = zeros_gm(n+1, n+1)

  	for ind1 in range(n+1):
    	for ind2 in range(n+1):
      		if ind1 <= ind2:
        		SGM[ind1, ind2] = inner_rjrk(ind1, ind2, lam)

  	SGM = symmetrize(SGM)

  	basis_mat = eye_gm(n+1)
  	o_basis_mat = zeros_gm(n+1, n+1)
  
  	o_basis_mat[0] = basis_mat[0]
  
  	print('Orthogonalizing Using Gram-Schmidt')
  	for r in tqdm.tqdm(range(1, n+1)):
      	u_r = basis_mat[r]
      	for i in range(r):
          	v_i = o_basis_mat[i]
          
          	proj = Polynomial.fast_inner(u_r, v_i, SGM)
          	norm = Polynomial.fast_inner(v_i, v_i, SGM)
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
                                          SGM)
              	o_basis_arr[i] = o_basis_mat[i]/gm.sqrt(norm)
          	return o_basis_arr
      	return np.array(o_basis_mat, dtype=np.float64)

  	return o_basis_mat


# num = 10
# frac = 1
# sops = generate_symm_ops(num, normalized = 1, frac = frac)
# print(sops)
# arr1 = sops[1]
# arr2 = sops[2]
# Polynomial.fast_inner(arr1, arr2, SGM)