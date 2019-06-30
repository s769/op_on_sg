from monomials import generate_T
from ops_main import generate_op
from Polynomial import Polynomial
from innerprods import lis2str
import scipy.io 
import numpy as np
#T = generate_T(7,31)
#print(T(1,1,0))
#scipy.io.savemat('../data/Tarray.mat', dict(T=T))
n = 5
k = 3
normalized = 1
lam=np.array([0])
ops_leg = generate_op(n, k, normalized, lam)
GM = Polynomial.GM[lis2str(lam)][:n+1, :n+1]
for i in range(20):
    for k in range(i+1):
        prod = Polynomial.fast_inner(ops_leg[i], ops_leg[k], GM)
        print("Inner product of "+str(i)+" "+str(k)+" "+str(prod))

