from Polynomial import *


def generate_op(n, k, normalized=True, lam=np.array([1])):
  Polynomial.build_condensed_GM(n+1, k, lam)
#   basis_mat = np.zeros((n+1, 3*n+3))
  
#   for i in range(n+1):
#     basis_mat[i, 3*i + k - 1] = 1
  basis_mat = np.eye(n+1)
  #o_basis_mat = np.zeros((n+1, 3*n+3))
  o_basis_mat = np.zeros((n+1, n+1))
  
  o_basis_mat[0] = basis_mat[0]
  for r in range(1, n+1):
    res = np.zeros(n+1)
    u_r = basis_mat[r]
    for i in range(r):
      v_i = o_basis_mat[i]
      le = len(u_r)
      GM = Polynomial.GM[lis2int(lam)][:le, :le]
      proj = Polynomial.fast_inner(u_r, v_i, GM)
      norm = Polynomial.fast_inner(v_i, v_i, GM)
      res += proj/norm*v_i
    o_basis_mat[r] = basis_mat[r] - res
    
  if normalized:
    for i in range(n+1):
      norm = Polynomial.fast_inner(o_basis_mat[i], o_basis_mat[i],
                                  GM)
      o_basis_mat[i] /= np.sqrt(norm)
  return o_basis_mat#, o_basis_mat[:, k-1::3]
  
j = 50
k = 3
normalized = 0
ops_sob = generate_op(j,k,normalized, lam=np.array([1]))
ops_leg = generate_op(j, k, normalized, lam=np.array([0]))
# #arr1 = ops[2]
# # arr2 = ops[3]
# # Polynomial.fast_inner(arr1, arr2, Polynomial.GM)]
# #print(ops)


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


scipy.io.savemat('..\\data\\coefs.mat', dict(ops=ops_sob))
# scipy.io.savemat('coefs2.mat', dict(ops=ops_leg))
# scipy.io.savemat('coefscomb.mat', dict(ops=ops))


# ops_sob

# ops = generate_op(j, k, 1, lam=np.ones(C))   
# ops
#ops_sob