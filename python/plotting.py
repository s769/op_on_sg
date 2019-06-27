import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from monomials import generate_T
import math

def fi(x, qi):
    '''
    This function is a contractive similarity of the plane centered
    at the point qi of dilation factor 1/2.

    support functions: none

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

    support functions: fi.m

    Args:
        x - 3^(m+1) by 2 matrix holding x and y coordinates
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
    from http://www.math.cornell.edu/~mhall/SGprograms/SG.m
    helper function for gaskplot

    '''

    q0 = np.array([np.cos(5*np.pi/12), np.sin(5*np.pi/12)])
    q1 = np.array([0, 0])
    q2 = np.array([np.cos(np.pi/12), np.sin(np.pi/12)])

    y = np.array([q0, q1, q2])

    for i in range(m):
        y = np.vstack((qcontract(y, q0), qcontract(y, q1), qcontract(y,q2)))
    
    
    ax = plt.axes(projection='3d')
    #ax.plot(y[:,0], y[:, 1], '.')


    return y, ax

def gaskplot(f, m):
    '''
    This function plots a function defined on the level m vertices of
    SG.

    support functions: SG.m

    Args:
        f - vector of size 3^(m+1) holding function values on vertices
        m - level of the function

    Returns:
        plots figure of the graph

    from http://www.math.cornell.edu/~mhall/SGprograms/gaskplot.m

    '''

    y, ax = SG(m)


    for k in range(3**m):
        xx = np.append(y[3*k:3*k+3, 0], y[3*k, 0])
        yy = np.append(y[3*k:3*k+3, 1], y[3*k, 1])
        zz = np.append(f[3*k:3*k+3], f[3*k])

        ax.plot(xx, yy, zz, 'b')

    plt.show()
    return

# m = 7

# T = generate_T(m, 1)
# f = T[2,:, -1]

# #f = np.ones(3**(m+1))

# gaskplot(f, m)

def getOmegas3(deg, level=7):
    T = generate_T(level, deg+1)
    W = T[2, :deg+2, :level+2]
    return W, T

def eval_antisymm(deg, level=7):
    W, T = getOmegas3(deg, level)
    coeff = W[deg]

    q = np.zeros(3**(level+1))

    for k in range(deg+1):
        if math.isclose(coeff[k], 0, abs_tol=1e-8):
            q += coeff[k]*T[2, :, k]

    return q

def plot_antisymm(num, level=7):

    for j in range(num):
        p = eval_antisymm(j, level)
        gaskplot(p, level)
        
    return

plot_antisymm(2)