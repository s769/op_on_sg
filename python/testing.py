import numpy as np
from monomials import *
from matplotlib import pyplot as plt
import gmpy2 as gm
import sympy as sp
from scipy.stats import linregress
import scipy.io
from plotting import plot_monomial, plot_op
from util import rotate_address, alternate_address
from symmetric import generate_symm_ops
from ops_main import generate_op_GS, generate_f
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
# def converttostr(arr):
#     string = ""
#     for j in range(len(arr)):
#         string = string + str(arr[j])
#     return string 

# def makerandompoints(level, degree):
    
#     addr = [""]
#     for i in range(3*(degree+1)):
#         while True:
#             flag = False
#             address = np.random.randint(3,size=level+1)
#             altaddress = alternate_address(level,address)
#             for k in range(len(addr)):
#                 if converttostr(address) == addr[k] or converttostr(altaddress) == addr[k]:
#                     flag = True
#             if flag == False:
#                 break
#         addr.append(converttostr(address))
#     addr = addr[1:len(addr)]
#     return addr

# def makerandommatrix(level, degree):
#     IMatrix = sp.zeros(3*(degree + 1))
#     addresses = makerandompoints(level, degree)
#     for i in range(3*degree + 3):
#         for j in range(3*degree + 3):
#             IMatrix[i,j] = p_jk(addresses[i], j//3 , (j%3)+1)
#     return IMatrix.det()

# counter = 0 
# for i in range(999):
#     if makerandommatrix(2,2) == 0:
#         counter = counter + 1

# print((counter/999))



# sob4_deg20 = generate_op_GS(20,3,False,np.array([1,1,1,1]),False)
# sob5_deg20 = generate_op_GS(20,3,False,np.array([1,1,1,1,1]),False)

# scipy.io.savemat('../data/sob4_deg20.mat', dict(coefs=sob4_deg20))
# scipy.io.savemat('../data/sob5_deg20.mat', dict(coefs=sob5_deg20))
# np.savez('../data/sob4_deg20', coefs=sob4_deg20)
# np.savez('../data/sob5_deg20', coefs=sob5_deg20)

f1_deg20 = generate_f(20,1,frac=False)
f2_deg20 = generate_f(20,2,frac=False)
f3_deg20 = generate_f(20,3,frac=False)


scipy.io.savemat('../data/f_1deg20.mat', dict(coefs=f1_deg20))
scipy.io.savemat('../data/f_2deg20.mat', dict(coefs=f2_deg20))
scipy.io.savemat('../data/f_3deg20.mat', dict(coefs=f3_deg20))

np.savez('../data/f_1deg20', coefs=f1_deg20)
np.savez('../data/f_2deg20', coefs=f2_deg20)
np.savez('../data/f_3deg20', coefs=f3_deg20)