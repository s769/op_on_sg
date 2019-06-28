from monomials import generate_T
import scipy.io 

T = generate_T(7,19)
print(T(1,1,0))
scipy.io.savemat('../data/Tarray.mat', dict(T=T))