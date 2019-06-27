from monomials import generate_T
import scipy.io 

T = generate_T(7,19)
scipy.io.savemat('../data/Tarray.mat', dict(T=T))