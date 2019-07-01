import numpy as np
import sympy as sp
from sympy import Rational as Rat

'''
This file computes recursively the alpha, beta, gamma, and eta coefficients 
from the Kasso, Tuley paper OP on SG.

The recursions are memoized (meaning they store the previously computed 
values in dictionaries) for increased performance. These functions are 
faster when you require specific values of the coefficients. 

For large computations that use all the coefficients for multiple times, 
it will be better to use the functions alpha_array, beta_array, etc..
'''

# This is a function that can be used to change the np.array float precision to 
#   float32 (double)


def array(*args, **kwargs):
    """
    This is a function that can be used to change the np.array float 
        precision to float32 (double)
    """
    kwargs.setdefault("dtype", np.float32)
    return np.array(*args, **kwargs)


def mem(func):
    '''
    This is a decorator function that memoizes the recursive functions. 
    This function is used as a decorator in multiple places thoughout the 
        code files.
    '''
    cache = dict()

    def memoized_func(*args):
        if args in cache:
            return cache[args]
        result = func(*args)
        cache[args] = result
        return result

    return memoized_func


'''
The follwing five functions compute single values of alpha, beta, gamma, eta, 
    and alpha prime from the Kasso, Tuley paper.
'''


@mem
def alpha(j):
    '''Calculates alpha_j

    Args: 
    j: index of coefficient

    Returns:
    alpha_j value
    '''
    if j == 0:
        return 1
    if j == 1:
        return Rat(1,6)
    res = 0
    for l in range(1, j):
        res += alpha(j-l)*alpha(l)
    return Rat(res*4,(5**j - 5))


@mem
def beta(j):
    '''Calculates beta_j

    Args: 
    j: index of coefficient

    Returns:
    beta_j value
    '''
    if j == 0:
        return Rat(-1,2)
    res = 0
    for l in range(j):
        res += (3*5**(j-l) - 5**(l+1) + 6)*alpha(j-l)*beta(l)
    return Rat(res*2,(15*(5**j - 1)))


@mem
def gamma(j):
    '''Calculates gamma_j

    Args: 
    j: index of coefficient

    Returns:
    gamma_j value
    '''
    return 3*alpha(j+1)


@mem
def eta(j):
    '''Calculates eta_j

    Args: 
    j: index of coefficient

    Returns:
    eta_j value
    '''
    if j == 0:
        return 0
    res = Rat(alpha(j)*(5**j + 1),2)
    for l in range(j):
        res += 2*eta(l)*beta(j-l)
    return res


@mem
def ap(j):
    '''Calculates alpha'_j

    Args: 
    j: index of coefficient

    Returns:
    alpha'_j value
    '''
    if j == 0:
        return Rat(1,2)
    return alpha(j)


# Here we calculate an array of values listed above.
# These are more efficient if we want to use all the values for
#   multiple times.

def alpha_array(max_order):
    """Calculates the array of alpha_j up to some order

    Args:
    max_order: Maximum order of j to be calculated (should be >= 1)

    Returns: 
    alpha_arr: np.array of length (max_order+1) containing the first 
        alpha values up to alpha_{max_order}.
    """
    alpha_arr = np.zeros((max_order + 1))
    alpha_arr[0] = 1
    alpha_arr[1] = 1 / 6

    for j in range(2, max_order+1):
        for l in range(1, j):
            alpha_arr[j] = alpha_arr[j] + alpha_arr[j-l] * alpha_arr[l]
        alpha_arr[j] = alpha_arr[j] * 4 / (5 ** j - 5)

    return alpha_arr


def beta_array(max_order):
    """Calculates the array of beta_j up to some order

    Args:
    max_order: Maximum order of j to be calculated (should be >= 1)

    Returns: 
    beta_arr: np.array of length (max_order+1) containing the first 
        alpha values up to beta_{max_order}.
    """

    # alpha's are used in the calculation of beta, so we first calculate
    #   these
    alpha_arr = alpha_array(max_order+2)

    # initialize values of beta_arr
    beta_arr = np.zeros((max_order + 1))
    beta_arr[0] = -1/2

    for j in range(1, max_order+1):
        for l in range(0, j):
            beta_arr[j] = beta_arr[j] + (3 * (5 ** (j-l)) -
                            5 ** (l + 1) + 6) * alpha_arr[j-l] * beta_arr[l]
        beta_arr[j] = beta_arr[j] * 2 / (15 * (5 ** j - 1))

    return beta_arr


def gamma_array(max_order):
    """Calculates the array of gamma_j up to some order

    Args:
    max_order: Maximum order of j to be calculated (should be >= 1)

    Returns: 
    gamma_arr: np.array of length (max_order+1) containing the first 
        alpha values up to gamma_{max_order}.
    """

    alpha_arr = alpha_array(max_order+2)
    gamma_arr = 3 * alpha_arr[1:max_order+2]
    return gamma_arr
