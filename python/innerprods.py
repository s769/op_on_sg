from recursions import *

import functools

# def inner_j1k1(j, k, energy=False):
#   if energy:
#     res = 2*alpha(k)*eta(j)
#     ms = min(j, k)
#     s1 = 0
#     for l in range(j-ms, j+1):
#       s1 += alpha(j-l)*eta(k+l+1) - alpha(k+l+1)*eta(j-l)
#     s2 = 0
#     for l in range(j-ms-1, j):
#       s2 += alpha(j-l-1)*eta(k+l+1) - alpha(k+l+1)*eta(j-l-1)

#     return res + 2*(s1-s2)
  

#   ms = min(j, k)
#   s1 = 0
#   for l in range(j-ms, j+1):
#     s1 += alpha(j-l)*eta(k+l+1) - alpha(k+l+1)*eta(j-l)
#   s2 = 0
#   for l in range(j-ms-1, j):
#     s2 += alpha(j-l-1)*n(k+l) - alpha(k+l)*n(j-l-1)
  
#   return 2*(s1+s2)


def inner0_j1k1(j, k):
  ms = min(j, k)
  s1 = 0
  for l in range(j-ms, j+1):
    s1 += alpha(j-l)*eta(k+l+1) - alpha(k+l+1)*eta(j-l)
  return 2*s1

# def inner_j2k2(j, k, energy=False):
#   if energy:
#     res = -2*beta(k)*alpha(j)
#     s1 = 0
#     ms = min(j, k)
#     for l in range(j-ms, j+1):
#       s1 += beta(j-l)*alpha(k+l+1) - beta(k+l+1)*alpha(j-l)
#     s2 = 0
#     for l in range(j-ms-1, j):
#       s2 += beta(j-l-1)*alpha(k+l+1) - beta(k+l+1)*alpha(j-l-1)

#     return res - 2*(s1-s2)
  
#   ms = min(j, k)
#   s1 = 0
#   for l in range(j-ms, j+1):
#     s1 += beta(j-l)*alpha(k+l+1) - beta(k+l+1)*alpha(j-l)
#   s2 = 0
#   for l in range(j-ms-1, j):
#     s2 += beta(j-l-1)*alpha(k+l) - beta(k+l)*alpha(j-l-1)
  
#   return -2*(s1+s2)

def inner0_j2k2(j, k):
  ms = min(j, k)
  s1 = 0
  for l in range(j-ms, j+1):
    s1 += beta(j-l)*alpha(k+l+1) - beta(k+l+1)*ap(j-l)
  return -2*s1

# def inner_j3k3(j, k, energy=False):
#   if energy:
#     res = 6*gamma(k)*eta(k)
#     ms = min(j, k)
#     s1 = 0
#     for l in range(j-ms, j+1):
#       s1 += alpha(j-l+1)*eta(k+l+2) - alpha(k+l+2)*eta(j-l+1)
#     s2 = 0
#     for l in range(j-ms-1, j):
#       s2 += alpha(j-l)*eta(k+l+2) - alpha(k+l+2)*eta(j-l)

#     return res + 18*(s1-s2)
  
#   ms = min(j, k)
#   s1 = 0
#   for l in range(j-ms, j+1):
#     s1 += alpha(j-l+1)*eta(k+l+2) - alpha(k+l+2)*eta(j-l+1)
#   s2 = 0
#   for l in range(j-ms-1, j):
#     s2 += alpha(j-l)*eta(k+l+1) - alpha(k+l+1)*eta(j-l)
  
#   return 18*(s1+s2)

def inner0_j3k3(j, k):
  ms = min(j, k)
  s1 = 0
  for l in range(j-ms, j+1):
    s1 += alpha(j-l+1)*eta(k+l+2) - alpha(k+l+2)*eta(j-l+1)
  return 18*s1

# def inner_j1k2(j, k, energy=False):
#   if energy:
#     res = 2*beta(k)*eta(k)
#     s1 = 0
#     for l in range(j+1):
#       s1 += alpha(j-l)*alpha(k+1+l) + beta(k+1+l)*eta(j-l)
#     s2 = 0  
#     for l in range(j):
#       s2 += alpha(j-l-1)*alpha(k+1+l) + beta(k+1+l)*eta(j-l-1)

#     return res - 2*(s1-s2)
  
#   s1 = 0
#   for l in range(j+1):
#     s1 += alpha(j-l)*alpha(k+l+1) + beta(k+l+1)*eta(j-l)
#   s2 = 0
#   for l in range(j):
#     s2 += alpha(j-l-1)*alpha(k+l) + beta(k+l)*eta(j-l-1)
  
#   return -2*(s1+s2)

def inner0_j1k2(j, k):
  s1 = 0
  for l in range(j+1):
    s1 += alpha(j-l)*alpha(k+l+1) + beta(k+l+1)*eta(j-l)
  return -2*s1
  
def symmetrize(arr):
  return arr + arr.T - np.diag(arr.diagonal())

def lis2int(lis):
  return functools.reduce(lambda total, d: 10 * total + d, lis, 0) 

inner0_j1k3 = lambda j, k: 0
inner0_j2k3 = lambda j, k: 0

inner_dict = {(1,1):inner0_j1k1, (2,2):inner0_j2k2, (3,3):inner0_j3k3,\
             (1,2):inner0_j1k2, (1,3):inner0_j1k3, (2,3): inner0_j2k3}


vals_dict = {1:alpha, 2:beta, 3:gamma}

dnpj2 = lambda j: -alpha(j)
dnpj3 = lambda j: 3*eta(j+1)

norm_dict = {1:eta, 2:dnpj2, 3:dnpj3}

