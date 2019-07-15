import numpy as np
import gmpy2 as gm


'''
This file computes recursively the alpha, beta, gamma, and eta 
    coefficients from the Kasso, Tuley paper OP on SG.

The recursions are memoized (meaning they store the previously computed 
    values in dictionaries) for increased performance. These functions 
    are faster when you require specific values of the coefficients. 

For large computations that use all the coefficients for multiple times, 
    it will be better to use the functions alpha_array, beta_array, etc.
'''


# This global variable can be used to change the recursion limit for 
#   better speed at high indices
j_max_rec = 1000


def zeros_gm(m, n):
    '''
    This function creates an m x n array of gmpy2 zeros. 
        It is used instead of np.zeros in places where rational 
        arithmetic is required.
    '''
    arr = np.array([[gm.mpz(0) for i in range(n)] for j in range(m)])
    return arr


def eye_gm(n):
    '''
    This function creates an n x n identity matrix using gmpy2 ones 
        and zeros. It is used instead of np.eye in places where rational 
        arithmetic is required.
    '''
    arr = zeros_gm(n, n)
    for i in range(n):
        arr[i, i] = gm.mpz(1)
    return arr

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

def mem2(func):
    '''
    This is a decorator function that memoizes recursive functions that 
        return arrays of values from index 0 to j. If the recursion has 
        already been called on a higher index, this decorator will cause 
        the function to automatically return those precomputed values.

    Example: Suppose mem2 decorates the function alpha_array(j), 
        which returns [alpha_0, ... , alpha_j]. Calling alpha_array(10) 
        will return [alpha_0, ... , alpha_10].
        Calling alpha_array(5) immediately after will return 
        [alpha_0, ... , alpha_10] without having to compute
        those values again.
    '''
    cache = dict()

    def memoized_func2(*args):
        check = False
        key_sp = None
        for key in cache.keys():
            if args <= key:
                check = True
                key_sp = key
                break
        if check:
            return cache[key_sp]
        result = func(*args)
        cache[args] = result
        return result
    return memoized_func2



# The follwing functions compute single values of alpha, beta, gamma, 
#     eta, and t from Calculus paper, and alpha prime from the Kasso, 
#     Tuley paper.



@mem
def alpha(j):
    '''Calculates alpha_j.

    Args: 
        j: index of coefficient

    Returns:
        value of alpha_j
    '''
    # For large index j, switch to array evaluation mode
    if j >= j_max_rec:
        return alpha_array(j)[j]
    
    # Otherwise, calculate directly based on (2.9) of Calculus I paper
    if j == 0:
        return 1
    if j == 1:
        return gm.mpq(1,6)
    res = 0
    for l in range(1, j):
        res += alpha(j-l)*alpha(l)
    return gm.mpq(res*4,(5**j - 5))


@mem
def beta(j):
    '''Calculates beta_j

    Args: 
        j: index of coefficient

    Returns:
        value of beta_j

    '''
    # For large index j, switch to array evaluation mode
    if j >= j_max_rec:
        return beta_array(j)[j]
    
    # Otherwise, calculate directly based on (2.11) of Calculus I paper
    if j == 0:
        return gm.mpq(-1,2)
    res = 0
    for l in range(j):
        res += (3*5**(j-l) - 5**(l+1) + 6)*alpha(j-l)*beta(l)
    return gm.mpq(res*2,(15*(5**j - 1)))


@mem
def gamma(j):
    '''Calculates gamma_j

    Args: 
    j: index of coefficient

    Returns:
    gamma_j value
    '''
    # For large index j, switch to array evaluation mode
    if j >= j_max_rec:
        return gamma_array(j)[j]
    
    # Otherwise, calculate directly based on (2.12) of Calculus I paper
    return 3*alpha(j+1)


@mem
def eta(j):
    '''Calculates eta_j

    Args: 
        j: index of coefficient

    Returns:
        value of eta_j
    '''
    # For large index j, switch to array evaluation mode
    if j >= j_max_rec:
        return eta_array(j)[j]
    
    # Otherwise, calculate directly based on (2.37) of Calculus I paper
    if j == 0:
        return 0
    res = gm.mpq(alpha(j)*(5**j + 1),2)
    for l in range(j):
        res += 2*eta(l)*beta(j-l)
    return res


@mem
def tau(j):
    '''Calculates tau_j

    Args: 
        j: index of coefficient

    Returns:
        value of tau_j
    '''
    # For large index j, switch to array evaluation mode
    if j >= j_max_rec:
        return tau_array(j)[j]
    
    # Otherwise, calculate directly based on (2.37) of Calculus I paper
    if j == 0:
        return gm.mpq(-1, 2)
    res = beta(j)
    for l in range(j):
        res -= 6*alpha(j+1-l)*tau(l)
    return res


@mem
def ap(j):
    '''Calculates alpha'_j

    Args: 
        j: index of coefficient

    Returns:
        alpha'_j value
    '''
    # For large index j, switch to array evaluation mode
    if j >= j_max_rec:
        return alpha_array(j)[j]
    
    # Otherwise, calculate directly based on (4) of Tuley paper
    if j == 0:
        return gm.mpq(1,2)
    return alpha(j)


# Here we calculate an array of values listed above.



@mem2
def alpha_array(max_order):
    """Calculates the array of alpha_j up to some order

    Args:
        max_order: Maximum order of j to be calculated (should be >= 1)

    Returns: 
        alpha_arr: np.array of length (max_order+1) containing the first 
            values of alpha up to alpha_{max_order}.
    """
    alpha_arr = zeros_gm(max_order + 1, 1)
    alpha_arr[0] = 1
    alpha_arr[1] = gm.mpq(1,6)

    for j in range(2, max_order+1):
        for l in range(1, j):
            alpha_arr[j] = alpha_arr[j] + alpha_arr[j-l] * alpha_arr[l]
        alpha_arr[j] = gm.mpq(alpha_arr[j] * 4 ,(5 ** j - 5))

    return alpha_arr


@mem2
def beta_array(max_order):
    """Calculates the array of beta_j up to some order

    Args:
        max_order: Maximum order of j to be calculated (should be >= 1)

    Returns: 
        beta_arr: np.array of length (max_order+1) containing the first 
            values of beta up to beta_{max_order}.
    """

    # alpha's are used in the calculation of beta, so we first calculate
    #   these
    alpha_arr = alpha_array(max_order+2)

    # initialize values of beta_arr
    beta_arr = zeros_gm(max_order + 1, 1)
    beta_arr[0] = gm.mpq(-1,2)

    for j in range(1, max_order+1):
        for l in range(j):
            beta_arr[j] = beta_arr[j] + (3 * (5 ** (j-l)) -
                            5 ** (l + 1) + 6) * alpha_arr[j-l] * beta_arr[l]
        beta_arr[j] = beta_arr[j] * gm.mpq(2 , (15 * (5 ** j - 1)))

    return beta_arr


@mem2
def gamma_array(max_order):
    """Calculates the array of gamma_j up to some order

    Args:
        max_order: Maximum order of j to be calculated (should be >= 1)

    Returns: 
        gamma_arr: np.array of length (max_order+1) containing the first 
            values of gamma up to gamma_{max_order}.
    """

    alpha_arr = alpha_array(max_order+2)
    gamma_arr = 3 * alpha_arr
    return np.array(gamma_arr[1:])


@mem2
def eta_array(max_order):
    """Calculates the array of eta_j up to some order

    Args:
        max_order: Maximum order of j to be calculated (should be >= 1)

    Returns: 
        eta_arr: np.array of length (max_order+1) containing the first 
            values of eta up to eta_{max_order}.
    """
    eta_arr = zeros_gm(max_order + 1, 1)
    eta_arr[0] = 0
    alpha_arr = alpha_array(max_order)
    beta_arr = beta_array(max_order)

    for j in range(1, max_order + 1):
        res = alpha_arr[j]*gm.mpq((5**j + 1), 2)
        for l in range(j):
            res += 2*eta_arr[l]*beta_arr[j-l]
        eta_arr[j] = res
    return eta_arr


@mem2
def tau_array(max_order):
    """Calculates the array of tau_j up to some order (this is 
        originally called t in the Calculus paper)

    Args:
        max_order: Maximum order of j to be calculated (should be >= 1)

    Returns: 
        tau_arr: np.array of length (max_order+1) containing the first 
            values of tau up to tau_{max_order}.
    """
    tau_arr = zeros_gm(max_order + 1, 1)
    tau_arr[0] = gm.mpq(-1, 2)

    alpha_arr = alpha_array(max_order+1)
    beta_arr = beta_array(max_order)

    for j in range(1, max_order + 1):
        res = beta_arr[j]
        for l in range(j):
            res -= 6*alpha_arr[j+1-l]*tau_arr[l]
        tau_arr[j] = res
    return tau_arr


@mem2
def ap_array(max_order):
    """Calculates the array of alpha'_j up to some order

    Args:
        max_order: Maximum order of j to be calculated (should be >= 1)

    Returns: 
        ap_arr: np.array of length (max_order+1) containing the first 
            values of alpha' up to alpha'_{max_order}.
    """
    ap_arr = alpha_array(max_order)
    ap_arr[0] = gm.mpq(1, 2)
    return ap_arr

