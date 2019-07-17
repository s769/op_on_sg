import numpy as np
from monomials import big_recursion, generate_T, generate_T_symmetric
from matplotlib import pyplot as plt
import gmpy2 as gm
from scipy.stats import linregress
import scipy.io
from plotting import plot_monomial, plot_op
from util import rotate_address
from symmetric import generate_symm_ops
# j = 50
# ar, br, pr, qr = big_recursion(j)
# a = np.array([float(x) for x in ar])
# b = np.array([float(x) for x in br])
# p = np.array([float(x) for x in pr])
# q = np.array([float(x) for x in qr])


# # plt.plot(np.log(np.abs(p)), 'bo', label='$p_j$')
# # plt.plot(np.log(np.abs(q)), 'rx', label='$q_j$')
# # plt.plot(np.log(np.abs(a)), 'g+', label='$a_j$')
# # plt.plot(np.log(np.abs(b)), 'k-', label='$b_j$')

# print(pr[1:]/pr[:-1] - qr[1:]/qr[:-1])
# # print()

# plt.plot(p[1:]/p[:-1], 'go')
# plt.plot(q[1:]/q[:-1], 'rx')

# plt.legend()

# indices = np.arange(j+1)

# # print(linregress(indices[1:], p[1:]))
# # print(linregress(indices[1:], q[1:]))


# plt.show()

# T = generate_T(7, 20, frac=False)

# scipy.io.savemat('../data/T20', dict(T=T))
# np.savez('../data/T20', T=T)

# ST = generate_T_symmetric(7, 20, frac=0, T=('../data/T20.npz', 'T'))
# #print(rotate_address(7, [0,0,0,0,0,0,0,0], 1))
# # scipy.io.savemat('../data/Tsymm20', dict(ST=ST))

# np.savez('../data/Tsymm20', ST=ST)

# # T = np.load('../data/T20.npz')['T']
# ST = np.load('../data/Tsymm20.npz')['ST']
# # # scipy.io.savemat('../data/T20.mat', dict(T=T))
# scipy.io.savemat('../data/Tsymm20.mat', dict(ST=ST))

# ops_sym = generate_symm_ops(20, frac=False)
# np.savez('../data/symops20', ops=ops_sym)

#plot_op(3, 3, T=('../data/Tsymm20.npz', 'ST'), coefs=('../data/symops20.npz', 'ops'), symm=True)

#plot_monomial(3, 3, T=('../data/Tsymm20.npz', 'ST'), symm=True)
