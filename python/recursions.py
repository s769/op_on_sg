from imports import *

def array(*args, **kwargs):
    kwargs.setdefault("dtype", np.float32)
    return np.array(*args, **kwargs)
  
def mem(func):
    cache = dict()

    def memoized_func(*args):
        if args in cache:
            return cache[args]
        result = func(*args)
        cache[args] = result
        return result

    return memoized_func

@mem
def alpha(j):
  if j==0: return 1
  if j==1: return 1/6
  res = 0
  for l in range(1, j):
    res += alpha(j-l)*alpha(l)
  return res*4/(5**j - 5)


@mem
def beta(j):
  if j==0: return -1/2
  res = 0
  for l in range(j):
    res += (3*5**(j-l) - 5**(l+1) + 6)*alpha(j-l)*beta(l)
  return res*2/(15*(5**j - 1))

@mem
def gamma(j):
  return 3*alpha(j+1)

@mem
def eta(j):
  if j==0: return 0
  res = alpha(j)*(5**j + 1)/2
  for l in range(j):
    res += 2*eta(l)*beta(j-l)
  return res

@mem
def ap(j):
  if j==0: return 1/2
  return alpha(j)