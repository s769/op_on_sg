from monomials import generate_T
from ops_main import generate_op
from Polynomial import Polynomial
from innerprods import lis2str
import scipy.io 
import numpy as np
import tqdm
from recursions import alpha, beta, gamma
T = generate_T(2, 2, frac=False)

#print(T)
#scipy.io.savemat('../data/Tarray.mat', dict(T=T))
# n = 5
# k = 3
# normalized = 1
# lam=np.array([0])
# ops_leg = generate_op(n, k, normalized, lam)
# GM = Polynomial.GM[lis2str(lam)][:n+1, :n+1]
# for i in range(20):
#     for k in range(i+1):
#         prod = Polynomial.fast_inner(ops_leg[i], ops_leg[k], GM)
#         print("Inner product of "+str(i)+" "+str(k)+" "+str(prod))

# def zeta(j, k, omega=None):
#   if omega is None:
#     omega = generate_op(j, 1, 0, lam=np.array([1]), frac=0)
  
#   if k == 1:
#     func = alpha
#   elif k == 2:
#     func = beta
#   elif k == 3:
#     func = gamma
#   res = 0
#   for l in range(j+1):
#     res += omega[j+1, l]*func(l+1)
    
#   return 2*res

# n = 200
# omega = generate_op(n, 1, 0, lam=np.array([0]), frac=0)
# norms = np.zeros(n+1)
# for i in range(len(omega)):
#   norms[i] = Polynomial.fast_inner(omega[i], omega[i],\
#                                   Polynomial.GM[0])
# norms = np.sqrt(norms)
# #print(norms)

# arr = np.zeros(n)

# for i in range(n):
#   arr[i] = zeta(i, 1, omega)/norms[i]
#   #arr[i] = zeta(i, 1, omega)
#   print(str(i)+':'+str(arr[i]))
# #print(np.sum(arr))
# for i in range(n-3):
#   num = (arr[i+3] - arr[i+2])/(arr[i+2]-arr[i+1])
#   den = (arr[i+2] - arr[i+1])/(arr[i+1]-arr[i])
#   #print(str(i)+':', arr[i+1]/arr[i])
#   #print(str(i)+':', np.log(np.abs(arr[i+1]/arr[i])))
#   #print(str(i) + ':', np.log(np.abs(num))/np.log(np.abs(den)))

# #print(arr)