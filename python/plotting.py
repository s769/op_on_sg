import math
import itertools

import numpy as np
import gmpy2 as gm
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

from monomials import generate_T, generate_W, generate_T_symmetric
from ops_main import generate_op
from symmetric import generate_symm_ops
from recursions import zeros_gm, eye_gm
import math
import itertools
import gmpy2 as gm
import tqdm


### GENERAL PLOTTING METHODS

def fi(x, qi):
    '''
    This function is a contractive similarity of the plane centered
    at the point qi of dilation factor 1/2.

    Args:
        x - point in the plane
        qi - point toward which to contract distances by 1/2

    Returns:
        evaluates the similarity

    from http://www.math.cornell.edu/~mhall/SGprograms/fi.m

    '''
    return qi + 0.5*(x-qi)

def qcontract(x, q):
    '''
    This function takes in a column of coordinate pairs and 
    computes their image under the similarity which contracts
    distance to the point q by a factor of 1/2.

    Args:
        x - (3^(m+1))*2 matrix holding x and y coordinates
        n - the number of points

    Returns: 
        y - coordinates of the n images

    from http://www.math.cornell.edu/~mhall/SGprograms/qcontract.m

    '''
    n = x.shape[0]
    y = np.zeros(x.shape)
    for i in range(n):
        y[i] = fi(x[i], q)
    return y

def SG(m):
    '''
    This function evaluates the coordinates of all points on the gasket
    up to a certain level V_m.

    Args:
        m - number of levels V_m we would like to find coordinates for
    
    Returns:
        y - ((3^m+3)/2)*2 matrix holding coordinate values

    from http://www.math.cornell.edu/~mhall/SGprograms/SG.m
    '''
    # Give coordinates for points in V_0
    q0 = np.array([np.cos(5*np.pi/12), np.sin(5*np.pi/12)])
    q1 = np.array([0, 0])
    q2 = np.array([np.cos(np.pi/12), np.sin(np.pi/12)])
    y = np.array([q0, q1, q2])
    
    # Evaluate all coordinate points 
    for _ in itertools.repeat(None, m):
        y = np.vstack((qcontract(y, q0), qcontract(y, q1), qcontract(y,q2)))
    
    # This plots the coordinates of points of V_m, if needed
    # ax.plot(y[:,0], y[:, 1], '.')
    return y

def gaskplot(f, m, ax):
    '''
    This function plots a function defined on the level m vertices of
    SG.

    Args:
        f - vector of size 3^(m+1) holding function values on vertices
        m - level of the function

    Returns:
        plots figure of the graph

    from http://www.math.cornell.edu/~mhall/SGprograms/gaskplot.m
    '''
    # Generate coordinates
    y = SG(m)

    for k in range(3**m):
        # Arrange the coordinates and the points
        xx = np.append(y[3*k:3*k+3, 0], y[3*k, 0])
        yy = np.append(y[3*k:3*k+3, 1], y[3*k, 1])
        zz = np.append(f[3*k:3*k+3], f[3*k])
        # Add the points to the plot
        ax.plot(xx, yy, zz, 'b')
    return


# METHODS FOR PLOTTING MONOMIALS

def plot_general(T, level=7):
    """
    Plot a general function given input

    Args: 
        num - Number of monomials we would like to plot
        k - Type of Monomial (k = 0, 1, 2)
        level - The level we would like to plot each monomial
    
    Returns: 
        figures of the SOP of type k, from P_{num-1, k} down to P_{0, k}.
    """
    plt.figure()
    ax = plt.axes(projection='3d')
    gaskplot(T, level, ax)
    plt.show()    
    return


def plot_monomial(num, k, level=7, T=None, symm=False):
    """
    Plot the Monomials

    Args: 
        num - Number of monomials we would like to plot
        k - Type of Monomial (k = 0, 1, 2)
        level - The level we would like to plot each monomial
    
    Returns: 
        figures of the SOP of type k, from P_{num-1, k} down to P_{0, k}.
    """
    if T is None:
        if symm:
            T = generate_T_symmetric(level, num, frac=False)
        else:
            T = generate_T(level, num, frac=False)
    elif isinstance(T, (tuple)):
        filename, arr = T
        T = np.load(filename, allow_pickle=False)[arr]

    for j in range(num):
        plt.figure()
        ax = plt.axes(projection='3d')
        p = T[:, j] if symm else T[k-1, :, j]
        gaskplot(p, level, ax)
    plt.show()    
    return


# METHODS FOR PLOTTING EASY BASIS

def plot_easy_basis(num, k, level=7, W=None):
    """
    Plot the Easy Basis

    Args: 
        num - Number of monomials we would like to plot
        k - Type of Monomial (k = 0, 1, 2)
        level - The level we would like to plot each monomial
    
    Returns: 
        figures of the SOP of type k, from P_{num-1, k} down to P_{0, k}.
    """

    if W is None:
        W = generate_W(level, num, frac=False)
    elif isinstance(W, (tuple)):
        filename, arr = W
        W = np.load(filename, allow_pickle=False)[arr]

    for j in range(num):
        plt.figure()
        ax = plt.axes(projection='3d')
        p = W[k, :, j]
        gaskplot(p, level, ax)
    plt.show()    
    return


# METHODS FOR PLOTTING SOBOLEV ORTHOGONAL POLYNOMIALS

def getOmegas(deg, k, frac=True, coefs=None, symm=False):
    """
    Function that generates the Sobolev Orthogonal Polynomials

    Args:
        deg - highest degree of the Sobolev Orthogonal Polynomial
            we would like to find
        k - Type of Sobolev Orthogonal Polynomial (k = 1,2,3)
    
    Returns:
        W - (deg+2)*(deg+2) matrix, representing the coefficients of 
            the Sobolev Orthogonal Polynomial of order 0 - deg+1
    """
    if coefs is None:
        if symm:
            W = generate_symm_ops(deg, normalized=1, frac=frac)
        else:
            W = generate_op(deg, k, 1, frac=frac)
    elif isinstance(coefs, tuple):
        filename, arr = coefs
        W = np.load(filename, allow_pickle=frac)[arr]

    else:
        W = coefs 

        
    # Generate the Sobolev orthogonal polynomials

    return W[:deg+2, :deg+2]

def eval_op(deg, k, level=7, T=None, frac=True, coefs=None, symm=False):
    """
    Function that evaluates the Sobolev Orthogonal Polynomials

    Args:
        deg - Degree of SOP s_{j} we would like to evaluate
        T - Coefficient Matrix of the Orthogonal Polynomials
        k - Type of SOP (k = 1, 2, 3)
        level - Number of levels we want to evaluate the SOP

    Returns:
        q - Values of the SOP of type k at some level
    """

    if T is None:
        if symm:
            T = generate_T_symmetric(level, deg, frac=frac)
        else:
            T = generate_T(level, deg, frac=frac)
    elif isinstance(T, (tuple)):
        filename, arr = T
        T = np.load(filename, allow_pickle=frac)[arr]



    # Fetch the particular coefficients of SOP 

    W = getOmegas(deg, k, frac=frac, coefs=coefs, symm=symm)

    coefs = W[:deg+1,:deg+1]
   

    Tarr = T[:, :deg+1] if symm else T[k-1, :, :deg+1]
    dtype = object if frac else np.float64
    q = np.empty((deg+1, Tarr.shape[0]), dtype=dtype)
    print('Evaluating Orthogonal Polynomials')
    for i in tqdm.tqdm(range(deg+1)):
        q[i] = np.sum(coefs[i]*Tarr, axis=1)
    return q

def plot_op(num, k, level=7, T=None, coefs=None, symm=False):
    """
    Plot the Sobolev Orthogonal Polynomials

    Args: 
        num - Number of Antisymmetric SOPs we would like to plot
        k - Type of SOP (k = 0, 1, 2)
        level - The level we would like to plot each SOP
    
    Returns: 
        figures of the SOP of type k, from s_{num-1} down to s_{0}.
    """
    p = eval_op(num, k, level=level, T=T, frac=False, coefs=coefs, symm=symm)
    for j in range(num):
        plt.figure()
        ax = plt.axes(projection='3d')
        gaskplot(p[j], level, ax)
    plt.show()    
    return



