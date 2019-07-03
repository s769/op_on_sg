from monomials import generate_T
from ops_main import generate_op
from Polynomial import Polynomial
from innerprods import lis2str
import scipy.io 
import numpy as np
import tqdm
from recursions import alpha, beta, gamma
import sympy as sp
from sympy import Rational as Rat
#T = generate_T(7, 19, frac=False)
#print("Now I will save everything!")
# # #print(T)
#scipy.io.savemat('../data/Tarray.mat', dict(T=T))
# n = 20
# k = 3
# normalized = 1
# lam=np.array([0])
# ops_leg = generate_op(n, k, normalized, lam, frac=False)
# # # print("Now I will save everything!")
# # print(ops_leg)
# scipy.io.savemat('../data/leg20coefs.mat', dict(ops=ops_leg))
# GM = Polynomial.GM[lis2str(lam)][:n+1, :n+1]
# for i in range(20):
#     for k in range(i+1):
#         prod = Polynomial.fast_inner(ops_leg[i], ops_leg[k], GM)
#         print("Inner product of "+str(i)+" "+str(k)+" "+str(prod))

def zeta(j, k, omega=None):
   if omega is None:
     omega = generate_op(j, 1, 0, lam=np.array([0]), frac=1)

   if k == 1:
     func = alpha
   elif k == 2:
     func = beta
   elif k == 3:
     func = gamma
   res = Rat(0,1)
   for l in range(j+1):
         res += omega[j, l]*func(l+1)
   return 2*res

n = 20
omega = generate_op(n, 1, 0, lam=np.array([0]), frac=1)
norms = sp.zeros(n+1,1)
for i in range(omega.rows):
         norms[i] = Polynomial.fast_inner(omega[i,:].T, omega[i,:].T,\
                                   Polynomial.GM['0'])[0]



# #print(norms)

arr = sp.zeros(n, 1)

for i in range(n):

  arr[i] = zeta(i, 1, omega)/sp.sqrt(norms[i])

  print(str(i)+':'+str(arr[i].evalf()))
#print(np.sum(arr))
for i in range(n-3):
  num = (arr[i+3] - arr[i+2])/(arr[i+2]-arr[i+1])
  den = (arr[i+2] - arr[i+1])/(arr[i+1]-arr[i])
  #print(str(i)+':', (arr[i+1]/arr[i]).evalf())
  #print(str(i)+':', np.log(np.abs(arr[i+1]/arr[i])))
  #print(str(i) + ':', np.log(np.abs(num))/np.log(np.abs(den)))

# #print(arr)