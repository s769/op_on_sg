import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from util import address_from_index
from monomials import generate_T
from plotting import gaskplot, plot_general
from recursions import alpha


# We first present functions specifically for H1

# Method for evaluating any function in H1 and returning all values
def eval_poly_h1(a, b, c, d, e, f, level = 7):
    """
    Function that evaluates a polynomial in H1. We can represent any 
    polynomial in H1 as the following linear combination:
    u = a * P_01 + b * P_02 + c * P_03 + d * P_11 + e * P_12 + f * P_13

    Args:
        a - coefficient of P_01
        b - coefficient of P_02
        c - coefficient of P_03
        d - coefficient of P_11
        e - coefficient of P_12
        f - coefficient of P_13
        level - The level we would like to evaluate 
    
    Returns:
        poly_temp - a np.array of length 3^(level+1), containing the 
            values of the polynomial.
    """
    T = generate_T(level, 2, frac=False)
    P01_temp = T[0, :, 0]
    P02_temp = T[1, :, 0]
    P03_temp = T[2, :, 0]
    P11_temp = T[0, :, 1]
    P12_temp = T[1, :, 1]
    P13_temp = T[2, :, 1]
    poly_temp = a * P01_temp + b * P02_temp + c * P03_temp \
        + d * P11_temp + e * P12_temp + f * P13_temp
    return poly_temp


# Method for plotting any function in H1 given coefficients
def plot_h1(a, b, c, d, e, f, level = 7):
    """
    Function that plots a polynomial in H1. This follows from the 
    evaluation function above.

    Args:
        a - coefficient of P_01
        b - coefficient of P_02
        c - coefficient of P_03
        d - coefficient of P_11
        e - coefficient of P_12
        f - coefficient of P_13
        level - The level we would like to plot
    
    Returns:
        A plot of the polynomial
    """
    poly_temp = eval_poly_h1(a, b, c, d, e, f, level)
    plot_general(poly_temp, level)


# Plot a family of polynomials based on one varying coefficient
def plot_h1_family(start, end, num_points, k, level = 7):
    """
    Function that plots the polynomial of one certain family in H1, with
    different coefficients of the P_0k term.

    Args:
        start - starting value of evaluation
        end - ending value of evaluation
        num_points - number of points we would like to evaluate
        k - The type of polynomial family (k=1,2,3)
        level - The level we would like to plot
    
    Returns:
        A plot of graphs of polynomials in the same picture.
    """
    T = generate_T(level, 2, frac=False)
    P01_temp = T[0, :, 0]
    P02_temp = T[1, :, 0]
    P03_temp = T[2, :, 0]
    P11_temp = T[0, :, 1]
    P12_temp = T[1, :, 1]
    P13_temp = T[2, :, 1]

    plt.figure()
    ax = plt.axes(projection='3d')
    color_list = ['red', 'orange', 'yellow', 'green', 'cyan', 'blue', \
        'purple', 'black', 'brown', 'gray']

    for i in range(num_points):
        t = np.linspace(start, end, num_points)[i]
        if (k == 1):
            poly_temp = t * P01_temp + P11_temp
        elif (k == 2):
            poly_temp = t * P02_temp + P12_temp
        else:
            poly_temp = t * P03_temp + P13_temp
        gaskplot(poly_temp, level, ax, color=color_list[i])

    plt.show()    
    return


# Methods for finding value and position of multiple maximal points 
#   for a single function in h1
def extremal_val_h1(a, b, c, d, e, f, num_points, flag, level = 7):
    """
    This function finds the points in a single polynomial in H^{1}, 
        with one of the following properties:
        (1) maximal value (2) minimal value (3) largest absolute value
        (4) minimum absolute value 
        and returns their value and address.

    Args:
        a - coefficient of P_01
        b - coefficient of P_02
        c - coefficient of P_03
        d - coefficient of P_11
        e - coefficient of P_12
        f - coefficient of P_13
        num_points - the number of maximal/minimal points we would like
            to check
        flag - the kind of extremum we're trying to find. 1 represent 
            maximal value, 2 represent minimal value, and 3 represent
            maximal of absolute value.
        level - The level we would like to evaluate 
    
    Returns:
        Prints the top values of the address and value of extremal points
    """

    T = generate_T(level, 2, frac=False)
    P01_temp = T[0, :, 0]
    P02_temp = T[1, :, 0]
    P03_temp = T[2, :, 0]
    P11_temp = T[0, :, 1]
    P12_temp = T[1, :, 1]
    P13_temp = T[2, :, 1]

    poly_temp = a * P01_temp + b * P02_temp + c * P03_temp +\
        d * P11_temp + e * P12_temp + f * P13_temp

    extremal_val(poly_temp, num_points, flag, level)


def max_h1_val_family(start, end, num_points, k, level = 7):
    """
    Function that prints a list of supremum values for a certain family,
    given different coefficients.

    Args:
        start - starting value of evaluation
        end - ending value of evaluation
        num_points - number of points we would like to evaluate
        k - The type of polynomial family (k=1,2,3)
        level - The level we would like to plot

    Returns:
        For each coefficient, prints out the address and value of 
            extremal points.
    """
    T = generate_T(level, 2, frac=False)
    P01_temp = T[0, :, 0]
    P02_temp = T[1, :, 0]
    P03_temp = T[2, :, 0]
    P11_temp = T[0, :, 1]
    P12_temp = T[1, :, 1]
    P13_temp = T[2, :, 1]

    index = np.linspace(start, end, num_points)

    for i in range(num_points):
        t = index[i]
        if (k == 1):
            poly_temp = t * P01_temp + P11_temp
        elif (k == 2):
            poly_temp = t * P02_temp + P12_temp
        else:
            poly_temp = t * P03_temp + P13_temp
        poly_temp = np.abs(poly_temp)
        max_val = np.max(poly_temp)
        max_pos = np.argmax(poly_temp)
        max_add = address_from_index(level, max_pos+1)

        print()
        print()
        print('Coefficient is ', t)
        print('Maximum value is', max_val)
        print('Maximum achieved at address ', max_add)


# We then provide functions for general H^{n}

def extremal_val(poly_temp, num_points, flag, level = 7):
    """
    This function finds the points in a single polynomial, 
        with one of the following properties:
        (1) maximal value (2) minimal value (3) largest absolute value 
        and returns their value and address.

    Args:
        poly_temp - evaluation of the polynomial to some level
        num_points - the number of maximal/minimal points
        flag - the kind of extremum we're trying to find. 1 represent 
            maximal value, 2 represent minimal value, and 3 represent
            maximal of absolute value.
        level - The level we would like to evaluate 
    
    Returns:
        Prints the top values of the address and value of extremal points
    """

    # If we need to find the sup-norm, take the absolute value
    if (flag == 3 or flag == 4):
        poly_temp = np.abs(poly_temp)
    
    # If we need to find the maximum, we flip the sequence
    if (flag == 1 or flag == 3):
        ind_arrange = np.flip(np.argsort(poly_temp))
    else: 
        # If we need to find the minimum, no need to flip
        ind_arrange = np.argsort(poly_temp)

    # Find the top points
    for i in range(num_points):
        temp_add = address_from_index(level, ind_arrange[i]+1)
        print("Address ", i)
        print(temp_add)
        print("Index", i)
        print(ind_arrange[i])
        print("Value ", i)
        print(poly_temp[ind_arrange[i]])


# A function to evaluate polynomials in H^{n}
def eval_poly_hn(coeff_arr, level=7):
    """
    Function that evaluates the polynomial on H^{n}.

    Args:
        coeff_arr - n*3 ndarray storing the coefficients of the polynomial
            coeff_arr[j, 0] represent coefficient of P_{j1}
            coeff_arr[j, 1] represent coefficient of P_{j2}
            coeff_arr[j, 2] represent coefficient of P_{j3}
        level - The level we would like to evaluate
    
    Returns:
        An array representing the values of the polynomial
    """
    # Fetch the highest power to evaluate
    max_j = coeff_arr.shape[0]

    # Evaluate all the polynomials
    T = generate_T(level, max_j+2, frac=False)

    poly_temp = np.zeros(T.shape[1])
    for j in range(max_j):
        for k in range(3):
            poly_temp = poly_temp + coeff_arr[j, k] * T[k, :, j]

    return poly_temp


# Method for plotting any function in Hn given coefficients
def plot_hn(coeff_arr, level = 7):
    """
    Function that plots a polynomial in H1. This follows from the 
    evaluation function above.

    Args:
        coeff_arr - n*3 ndarray storing the coefficients of the polynomial
            coeff_arr[j, 0] represent coefficient of P_{j1}
            coeff_arr[j, 1] represent coefficient of P_{j2}
            coeff_arr[j, 2] represent coefficient of P_{j3}
        level - The level we would like to plot
    
    Returns:
        A plot of the polynomial
    """
    poly_temp = eval_poly_hn(coeff_arr, level)
    plot_general(poly_temp, level)


# Methods for finding value and position of multiple maximal points 
#   for a single function in H^{n}
def extremal_val_hn(coeff_arr, num_points, flag, level = 7):
    """
    This function finds the points in a single polynomial in H^{1}, 
        with one of the following properties:
        (1) maximal value (2) minimal value (3) largest absolute value
        (4) minimum absolute value 
        and returns their value and address.

    Args:
        coeff_arr - n*3 ndarray storing the coefficients of the polynomial
            coeff_arr[j, 0] represent coefficient of P_{j1}
            coeff_arr[j, 1] represent coefficient of P_{j2}
            coeff_arr[j, 2] represent coefficient of P_{j3}
        num_points - the number of maximal/minimal points we would like
            to evaluate
        flag - the kind of extremum we're trying to find. 1 represent 
            maximal value, 2 represent minimal value, 3 represent
            maximal of absolute value, 4 represent minimal of absolute
            value.
        level - The level we would like to evaluate 
    
    Returns:
        Prints the top values of the address and value of extremal points
    """

    poly_temp = eval_poly_hn(coeff_arr, level)
    extremal_val(poly_temp, num_points, flag, level)


def max_h2_val_family(start0, end0, num0, start1, end1, num1, num_max, \
    k, level=7, verbose_flag=False):
    """
    Function that prints a list of supremum values for a certain family,
    given different coefficients. [This essentially acts like a grid
    search of coefficients.]

    Args:
        start0 - starting value for P_0k coefficient
        end0 - ending value of P_0k coefficient
        num0 - number of points to evaluate for P_0k coefficient
        start1 - starting value for P_1k coefficient
        end1 - ending value of P_1k coefficient
        num1 - number of points to evaluate for P_1k coefficient
        num_max - maximum number of output coefficients with smallest 
            norm that we would like to see.
        k - The type of polynomial family (k=1,2,3)
        level - The level we would like to plot
        verbose - A flag that determines if we print best stats and 
            address during the iterations

    Returns:
        For each coefficient, prints out the address and value of 
            extremal points.
    """
    # Generate values of the polynomial
    T = generate_T(level, 3, frac=False)
    P01_temp = T[0, :, 0]
    P02_temp = T[1, :, 0]
    P03_temp = T[2, :, 0]
    P11_temp = T[0, :, 1]
    P12_temp = T[1, :, 1]
    P13_temp = T[2, :, 1]
    P21_temp = T[0, :, 2]
    P22_temp = T[1, :, 2]
    P23_temp = T[2, :, 2]

    # Generate space of coefficients we would like to grid search on
    index_0 = np.linspace(start0, end0, num0)
    index_1 = np.linspace(start1, end1, num1)

    # Storage of the values obtained at the best coefficients
    xs = np.zeros((num0, num1))
    ys = np.zeros((num0, num1))
    zs = np.zeros((num0, num1))
    pos_arr = np.zeros((num0, num1))

    # Main loop over all coefficients
    for i in range(num0):
        for j in range(num1):
            # Pick the coefficient value
            t0 = index_0[i]
            t1 = index_1[j]

            # Evaluate the polynomial
            if (k == 1):
                poly_temp = t0 * P01_temp + t1 * P11_temp + P21_temp
            elif (k == 2):
                poly_temp = t0 * P02_temp + t1 * P12_temp + P22_temp
            else:
                poly_temp = t0 * P03_temp + t1 * P13_temp + P23_temp
            
            # Find the max-norm of the polynomial
            poly_temp = np.abs(poly_temp)
            max_val = np.max(poly_temp)
            max_pos = np.argmax(poly_temp)
            max_add = address_from_index(level, max_pos+1)
            
            # Store the max-norm and the address
            xs[i, j] = t0
            ys[i, j] = t1
            zs[i, j] = max_val
            pos_arr[i, j] = max_pos
            
            # Print max-norm and address, if requested
            if verbose_flag:
                print()
                print()
                print('Coefficient of P0k is ', t0)
                print('Coefficient of P1k is ', t1)
                print('Maximum value is', max_val)
                print('Maximum achieved at address ', max_add)
    
    # Reshape the arrays to have the right shape
    xs = xs.reshape(num0*num1)
    ys = ys.reshape(num0*num1)
    zs = zs.reshape(num0*num1)

    # Sort the array based on the minimum max_norm
    index = np.argsort(zs)

    # Find the best coefficients and return their values
    for i in range(num_max):
        temp_ind = index[i]
        
        print("Coefficient of P0k with ranking ", i)
        print(xs[temp_ind])
        print("Coefficient of P1k with ranking ", i)
        print(ys[temp_ind])
        print("Value with ranking", i)
        print(zs[temp_ind])


# def max_hn_val_family(n, range_list, num_points, k, level = 7):
#     """
#     Function that prints a list of supremum values for a certain family,
#     given different coefficients.

#     Args:
#         n - The dimension of space H^{n} we're working with
#         range_list - A 2*n ndarray, with the i-th column representing
#             the lower and upper bound we search for the best value
#         num_points - A n dim array, each entry represents the number of 
#             points we would like to evaluate in a certain dimension
#         k - The type of polynomial family (k=1,2,3)
#         level - The level we would like to plot

#     Returns:
#         For each coefficient, prints out the address and value of 
#             extremal points.
#     """
    
#     T = generate_T(level, 3, frac=False)
#     P01_temp = T[0, :, 0]
#     P02_temp = T[1, :, 0]
#     P03_temp = T[2, :, 0]
#     P11_temp = T[0, :, 1]
#     P12_temp = T[1, :, 1]
#     P13_temp = T[2, :, 1]
#     P21_temp = T[0, :, 2]
#     P22_temp = T[1, :, 2]
#     P23_temp = T[2, :, 2]

#     index = np.linspace(start, end, num_points)

#     for i in range(num_points):
#         t = index[i]
#         if (k == 1):
#             poly_temp = t * P01_temp + P11_temp
#         elif (k == 2):
#             poly_temp = t * P02_temp + P12_temp
#         else:
#             poly_temp = t * P03_temp + P13_temp
#         poly_temp = np.abs(poly_temp)
#         max_val = np.max(poly_temp)
#         max_pos = np.argmax(poly_temp)
#         max_add = address_from_index(level, max_pos+1)

#         print()
#         print()
#         print('Coefficient is ', t)
#         print('Maximum value is', max_val)
#         print('Maximum achieved at address ', max_add)


#     num_points_total = np.prod(num_points)
#     for i in range(num_points_total):
#         # Fetch the index of the 

    

    
#     for i in range(n):
#         for j in range(num_points[i]):
#             if (k == 1):
#                 T[1,]



#     for i in range(num_points):
        
#         t = index[i]
#         if (k == 1):
#             poly_temp = t * P01_temp + P11_temp
#         elif (k == 2):
#             poly_temp = t * P02_temp + P12_temp
#         else:
#             poly_temp = t * P03_temp + P13_temp
#         poly_temp = np.abs(poly_temp)
#         max_val = np.max(poly_temp)
#         max_pos = np.argmax(poly_temp)
#         max_add = address_from_index(level, max_pos+1)

#         print()
#         print()
#         print('Coefficient is ', t)
#         print('Maximum value is', max_val)
#         print('Maximum achieved at address ', max_add)