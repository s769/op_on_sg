import numpy as np
from recursions import alpha, beta, gamma, eta, zeros, eye
from innerprods import inner0_j1k1, symmetrize
from Polynomial import Polynomial
import sympy as sp
from sympy import Rational as Rat
import tqdm
import gmpy2 as gm

def pnj1(nn, j, i):
  return int(j == 0) if i == nn else alpha(j)

def d_pnj1(nn, j, i):
  return 0 if i == nn else eta(j)


def inner_pnj1pmk1(nn, j, m, k):
  if nn == m: return inner0_j1k1(j,k)
  res = 0
  for l in range(j+1):
    s1 = 0
    for i in range(3):
      s1 += pnj1(nn, j-l, i)*d_pnj1(m, k+l+1, i)\
        -pnj1(m, k+l+1, i)*d_pnj1(nn, j-l, i)
    res += s1
  return res


def inner_rjrk(j,k, lam = 1):
  res = 0
  for nn in range(3):
    for m in range(3):
      res += inner_pnj1pmk1(nn, j, m, k) \
      + lam*inner_pnj1pmk1(nn, j-1, m, k-1)
  return res








def generate_symm_ops(n, normalized=False, frac=True):
  
  print('Building Gram Matrix ... this may take some time')
  SGM = zeros(n+1, n+1)

  for ind1 in range(n+1):
    for ind2 in range(n+1):
      if ind1 <= ind2:
        SGM[ind1, ind2] = inner_rjrk(ind1, ind2)

  SGM = symmetrize(SGM)
#   basis_mat = np.zeros((n+1, 3*n+3))

#   for i in range(n+1):
#     basis_mat[i, 3*i + k - 1] = 1
  basis_mat = eye(n+1)
  #o_basis_mat = np.zeros((n+1, 3*n+3))
  o_basis_mat = zeros(n+1, n+1)
  

  o_basis_mat[0,:] = basis_mat[0,:]
  print('Orthogonalizing Using Gram-Schmidt')
  for r in tqdm.tqdm(range(1, n+1)):
      #res = np.zeros(n+1)

      u_r = basis_mat[r,:]
      for i in range(r):
          v_i = o_basis_mat[i,:]
          proj = Polynomial.fast_inner(u_r.T, v_i.T, SGM)
          norm = Polynomial.fast_inner(v_i.T, v_i.T, SGM)
          u_r -= (proj/norm)*v_i
      o_basis_mat[r,:] = u_r#basis_mat[r] - res

  if normalized:
      print('Normalizing')
      for i in tqdm.tqdm(range(n+1)):
          norm = Polynomial.fast_inner(o_basis_mat[i,:].T, o_basis_mat[i,:].T,
                                        SGM)
          o_basis_mat[i,:] = o_basis_mat[i,:]/gm.sqrt(norm[0])

  return (o_basis_mat , SGM) if frac else (np.array(o_basis_mat).astype(np.float64), SGM)# , o_basis_mat[:, k-1::3]

num = 10
frac = 0
sops = generate_symm_ops(num, 1, frac)[0]
print(sops)
# arr1 = sops[1]
# arr2 = sops[2]
# Polynomial.fast_inner(arr1, arr2, SGM)