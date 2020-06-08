import numpy as np
from monomials import *
from matplotlib import pyplot as plt
import gmpy2 as gm
import sympy as sp
from scipy.stats import linregress
import scipy.io
from plotting import plot_monomial, plot_op, eval_op, gaskplot
from util import rotate_address, alternate_address
from symmetric import generate_symm_ops
from ops_main import generate_op_GS, generate_f, sob_ops_recursion
from recursions import gamma_array, alpha_array, beta_array
from Polynomial import Polynomial
from innerprods import lis2str
from scipy.stats import linregress
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

# sob2_deg20 = generate_op_GS(20,3,False,np.array([1,1]),False)
# sob3_deg20 = generate_op_GS(20,3,False,np.array([1,1,1]),False)
# sob4_deg20 = generate_op_GS(20,3,False,np.array([1,1,1,1]),False)
# sob5_deg20 = generate_op_GS(20,3,False,np.array([1,1,1,1,1]),False)

# scipy.io.savemat('../data/sob2_deg20.mat', dict(coefs=sob2_deg20))
# scipy.io.savemat('../data/sob3_deg20.mat', dict(coefs=sob3_deg20))
# scipy.io.savemat('../data/sob4_deg20.mat', dict(coefs=sob4_deg20))
# scipy.io.savemat('../data/sob5_deg20.mat', dict(coefs=sob5_deg20))
# np.savez('../data/sob2_deg20', coefs=sob2_deg20)
# np.savez('../data/sob3_deg20', coefs=sob3_deg20)
# np.savez('../data/sob4_deg20', coefs=sob4_deg20)
# np.savez('../data/sob5_deg20', coefs=sob5_deg20)

# f1_deg20 = generate_f(20,1,frac=False)
# f2_deg20 = generate_f(20,2,frac=False)
# f3_deg20 = generate_f(20,3)


# scipy.io.savemat('../data/f_1deg20.mat', dict(coefs=f1_deg20))
# scipy.io.savemat('../data/f_2deg20.mat', dict(coefs=f2_deg20))
# scipy.io.savemat('../data/f_3deg20.mat', dict(coefs=f3_deg20))

# np.save('../data/f_1deg20', f1_deg20)
# np.save('../data/f_2deg20', f2_deg20)
# np.save('../data/f_3deg20', f3_deg20)


# p11 = generate_op_GS(10,1)
# p12 = sob_ops_recursion(10,1)

# print(np.sum(np.abs(p11-p12)))

# p21 = generate_op_GS(20,1)
# p22 = sob_ops_recursion(20,1)

# print(np.sum(np.abs(p21-p22)))

# p31 = generate_op_GS(20,1)
# p32 = sob_ops_recursion(20,1)

# print(np.sum(np.abs(p31-p32)))

##### TESTING f RECURSION #####

# sob_coefs = generate_op_GS(20,3)
# f3_deg20 = generate_f(20,3)

# frac = True
# leg_omegas = generate_op_GS(20,3, lam=np.array([0]))
# normalized = False
# j = 20
# k = 3
    
# o_basis_mat = np.empty((j+1, j+1), dtype=object)
# f_mat = np.empty((j+1, j+1), dtype=object)
# print('Using Gram-Schmidt to generate initial Sobolev Polynomials')

# first_mat = generate_op_GS(1, k, normalized=False, frac=frac)
# const = gm.mpz(0) if frac else 0
# first_mat = np.pad(first_mat, ((0,0), (0, j-1)),'constant', constant_values=(const,))
# o_basis_mat[:2] = first_mat


# if k == 3: func_array = gamma_array
# if k == 2: func_array = beta_array
# if k == 1: func_array = alpha_array

# print('Generating values for f_n')
# func_arr = func_array(j+2)
# print('Building Gram Matrix for inner product caluclation.')
# Polynomial.build_condensed_GM(j+1, k, np.array([1]))
# GM = Polynomial.GM[lis2str(np.array([1]))][:j+1, :j+1]



# print('Using recursion to generate the rest of the Sobolev Polynomials')

# for ind in tqdm.tqdm(range(1,j), file=sys.stdout):
#     func_vec = func_arr[1:ind+2]
#     omega_vec = leg_omegas[ind, :ind+1]
#     zeta_ind = gm.mpq(-1,func_arr[0])*func_vec.dot(omega_vec)
#     f_ind = np.insert(omega_vec, 0, zeta_ind)
#     a_ind = Polynomial.fast_inner(f_ind, o_basis_mat[ind,:ind+2], GM[:ind+2, :ind+2])
#     b_ind = Polynomial.fast_inner(f_ind, o_basis_mat[ind-1,:ind+2], GM[:ind+2, :ind+2])
#     a_ind = gm.mpq(a_ind, Polynomial.fast_inner(o_basis_mat[ind, :ind+1], o_basis_mat[ind,:ind+1], GM[:ind+1, :ind+1]))
#     b_ind = gm.mpq(b_ind, Polynomial.fast_inner(o_basis_mat[ind-1, :ind], o_basis_mat[ind-1,:ind], GM[:ind, :ind]))
#     new_vec = f_ind - a_ind*o_basis_mat[ind, :ind+2] - b_ind*o_basis_mat[ind-1, :ind+2]

#     o_basis_mat[ind+1] = np.pad(new_vec, (0, j-ind-1), 'constant', constant_values=(const,))
#     f_mat[ind+1] = sob_coefs[ind+1] + a_ind*sob_coefs[ind] + b_ind*sob_coefs[ind-1]

# print(f3_deg20[:3])
# print(f_mat[2:5])
# print(np.sum(np.abs(f_mat[2:] - f3_deg20[2:])))

# lam_arr = ["../data/sob_lambda10.npz", "../data/sob_lambda100.npz", "../data/sob_lambda1000.npz", "../data/sob_lambda10000.npz"]
# f2_arr = eval_op(20,2,T=("../data/T20.npz", "T"), coefs=('../data/f_2deg20.npz', "coefs"))
# l_inf_arr = []
# j = 10
# for path in lam_arr:
#     s_arr = eval_op(20,2,T=("../data/T20.npz", "T"),coefs=(path, "coefs"))
#     l_inf_arr.append(np.max(np.abs(f2_arr[j-1]-s_arr[j])))
#     print(type((f2_arr[j-1]-s_arr[j])[0]))



# plt.plot([1,2,3,4], np.log10(l_inf_arr), "bo-")
# # slope, intercept, r_value, p_value, std_err = linregress(range(2,20), np.log10(l_inf_arr))
# # print("k = 2 slope: %f    intercept: %f" % (slope, intercept))
# # plt.title("k = 2 slope: %f    intercept: %f" % (slope, intercept))
# plt.show()


deg = 3


lam_array = [1000,10000,100000,1000000,10000000,100000000,1000000000,10000000000]
f3_arr = np.array(generate_f(20,3),dtype=np.float64)
f3_deg20 = generate_f(20,3,frac=False)


# # coeff_arr = []
# # vals_arr = []
f3_vals = eval_op(deg+2,3,T=("../data/T20.npz", "T"), coefs=f3_arr)
scipy.io.savemat('../data/f3new.mat', dict(coefs=f3_arr))
# # f3_load = np.load("../data/f_3deg20.npy")

# # print(np.sum(np.abs(f3_arr-f3_load)))
print(np.sum(np.abs(f3_arr-f3_deg20)))
# for lamb in lam_array:

#     sob_arr = generate_op_GS(deg+2,3,lam=np.array([lamb]))

#     # l_inf_arr = np.array([np.max(np.abs(f3_arr[i]-sob_arr[i])) for i in range(3,20)], dtype=np.float64)
#     # print("coeff lambda:", lamb)
#     # print(l_inf_arr)
#     # coeff_arr.append(np.max(np.abs(f3_arr[deg]-sob_arr[deg])))


#     sob_vals = np.array(eval_op(deg+2,3,T=("../data/T20.npz", "T"), coefs=sob_arr),dtype=np.float64)

#     plt.figure()
#     ax = plt.axes(projection='3d')
#     gaskplot(sob_vals[deg], 7, ax)
#     plt.title('lambda: {l}'.format(l=lamb))
    
    # l_inf_vals = np.array([np.max(np.abs(f3_vals[i]-sob_vals[i])) for i in range(3,20)], dtype=np.float64)
    # print("vals lambda:", lamb)
    # print(l_inf_vals)
    # vals_arr.append(np.max(np.abs(f3_vals[deg]-sob_vals[deg])))


plt.figure()
ax = plt.axes(projection='3d')
gaskplot(f3_vals[deg], 7, ax)
plt.title('$f_{3}$')

# plt.show()
# coeff_arr = np.log10(np.array(coeff_arr, dtype=np.float64))
# vals_arr = np.log10(np.array(vals_arr, dtype=np.float64))


# plt.figure()
# plt.plot(range(3,11), coeff_arr, "o-")
# slope, intercept, r_value, p_value, std_err = linregress(range(3,11), coeff_arr)
# plt.title("coefs slope: {sl} intercept: {int}".format(sl=slope, int=intercept))
# plt.figure()
# plt.plot(range(3,11), vals_arr, "o-")
# slope, intercept, r_value, p_value, std_err = linregress(range(3,11), vals_arr)
# plt.title('vals slope: {sl} intercept: {int}'.format(sl=slope, int=intercept))
plt.show()