import numpy as np
from matplotlib import pyplot as plt

from util import address_from_index
from monomials import generate_T
from plotting import gaskplot, plot_general
from recursions import alpha


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
    if (flag == 3):
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


# Methods for finding value and position of multiple maximal points 
#   for a single function in h1
def extremal_val_h1(num_points, a, b, c, d, e, f, flag, level = 7):
    """
    This function finds the points in a single polynomial in H^{1}, 
        with one of the following properties:
        (1) maximal value (2) minimal value (3) largest absolute value 
        and returns their value and address.

    Args:
        num_points - the number of maximal/minimal points
        a - coefficient of P_01
        b - coefficient of P_02
        c - coefficient of P_03
        d - coefficient of P_11
        e - coefficient of P_12
        f - coefficient of P_13
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


def poly_sup_division(j1, k1, j2, k2, ):
    """
    This function calculates the 
    """
    pass