import numpy as np
from monomials import big_recursion
from matplotlib import pyplot as plt
import gmpy2 as gm
from scipy.stats import linregress

j = 50
ar, br, pr, qr = big_recursion(j)
a = np.array([float(x) for x in ar])
b = np.array([float(x) for x in br])
p = np.array([float(x) for x in pr])
q = np.array([float(x) for x in qr])


# plt.plot(np.log(np.abs(p)), 'bo', label='$p_j$')
# plt.plot(np.log(np.abs(q)), 'rx', label='$q_j$')
# plt.plot(np.log(np.abs(a)), 'g+', label='$a_j$')
# plt.plot(np.log(np.abs(b)), 'k-', label='$b_j$')

print(pr[1:]/pr[:-1] - qr[1:]/qr[:-1])
# print()

plt.plot(p[1:]/p[:-1], 'go')
plt.plot(q[1:]/q[:-1], 'rx')

plt.legend()

indices = np.arange(j+1)

# print(linregress(indices[1:], p[1:]))
# print(linregress(indices[1:], q[1:]))


plt.show()
