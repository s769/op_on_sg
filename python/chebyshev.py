import numpy as np
from matplotlib import pyplot as plt

from util import address_from_index
from monomials import generate_T
from plotting import gaskplot, plot_general
from recursions import alpha


# Method for plotting any function in H0


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
        poly_temp - a np.array of length 3^(level+1), containing the 
            values of the polynomial.
    """
    poly_temp = eval_poly_h1(a, b, c, d, e, f, level)
    plot_general(poly_temp, level)


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
    color_list = ['red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'purple', 'black', 'brown', 'gray']

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


# Methods for finding numerical values of maximal/minimal points
def max_h1_family(start, end, num_points, k, level = 7):
    T = generate_T(level, 2, frac=False)
    P01_temp = T[0, :, 0]
    P02_temp = T[1, :, 0]
    P03_temp = T[2, :, 0]
    P11_temp = T[0, :, 1]
    P12_temp = T[1, :, 1]
    P13_temp = T[2, :, 1]

    index = np.linspace(start, end, num_points)
    max_list = np.zeros(num_points)

    for i in range(num_points):
        t = index[i]
        if (k == 1):
            poly_temp = t * P01_temp + P11_temp
        elif (k == 2):
            poly_temp = t * P02_temp + P12_temp
        else:
            poly_temp = t * P03_temp + P13_temp
        poly_temp = np.abs(poly_temp)
        max_list[i] = max(poly_temp)

    print()
    print()
    print(index)
    print(max_list)
    

# Methods for finding value and position of maximal points at a single
def max_h1_val(t, k, num_max, level = 7):
    # t is the constant denoting P_{1k}+tP_{0k}
    T = generate_T(level, 2, frac=False)
    P01_temp = T[0, :, 0]
    P02_temp = T[1, :, 0]
    P03_temp = T[2, :, 0]
    P11_temp = T[0, :, 1]
    P12_temp = T[1, :, 1]
    P13_temp = T[2, :, 1]

    if (k == 1):
        poly_temp = t * P01_temp + P11_temp
    elif (k == 2):
        poly_temp = t * P02_temp + P12_temp
    else:
        poly_temp = t * P03_temp + P13_temp
    poly_temp = np.abs(poly_temp)
    ind_arrange_1 = np.flip(np.argsort(poly_temp))

    for i in range(num_max):       
        temp_add = address_from_index(level, ind_arrange_1[i]+1)
        print("Address ", i)
        print(temp_add)
        print("Index", i)
        print(ind_arrange_1[i])
        print("Value ", i)
        print(poly_temp[ind_arrange_1[i]])





# Methods for finding value and position of maximal points at a single
def max_h1_val_family(start, end, num_points, k, level = 7):
    # t is the constant denoting P_{1k}+tP_{0k}
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