import numpy as np

from recursions import alpha_array, gamma_array

import sys
import gmpy2 as gm

# This file will implement helper functions written in the matlab files.

def address_from_index(level, index):
    """Computes address vector from the index.

    Returns the address vector of a point with certain TLR index in 
        a graph of SG at some fixed level.
    
    For point F_{w1} F_{w2} ... F_{wm} q_{k} (0<=wi<=2, 0<=k<=2), we can 
    represent it in two ways:
    1. Address Vector: [wm, ..., w1, k] (Note that this represents
        picking transform F_{wm}, ..., F_{w1} in this sequence and then 
        finding vertex q_{k} in the final triangle. This order is weird 
        but will prove to be useful)
    2. TLR Index: k = wm * 3^m + ... + w1 * 3 + k + 1

    This function provides transform 2 -> 1.

    Args:
        level: A nonnegative integer representing the level of SG we're 
            working with.
        index: A number in the range [1, 3^(level+1)] representing the 
            TLR index of some point in level m.
    
    Returns:
        address: np.array of length (level+1) representing the address 
            of the point specified by index in a particular level of SG.
    
    Example:
        address_from_index(1, 4) = [1, 0]
        address_from_index(2, 21) = [2, 1, 2]
    """

    # Initialize the vector to store the address
    v = np.zeros(level+1)

    # First decide the last term of the vector
    index = index - 1
    v[level] = index % 3
    index = int(index / 3)

    # Perform transformation similar to base 10 -> base 3
    for i in range(level):
        v[level - 1 - i] = index % 3
        index = int(index / 3)
    
    return v


def index_from_address(level, address):
    """Computes TLR Index from Address vector.

    Returns the TLR index of a point with certain address vector in 
        a graph of SG at some fixed level.
    
    For point F_{w1} F_{w2} ... F_{wm} q_{k} (0<=wi<=2, 0<=k<=2), we 
        can represent it in two ways:
    1. Address Vector: [wm, ..., w1, k] (Note that this represents
        picking transform F_{wm}, ..., F_{w1} in this sequence and then 
        finding vertex q_{l} in the final triangle)
    2. TLR Index: k = wm * 3^m + ... + w1 * 3 + k + 1

    This function provides transform 1 -> 2.

    Args:
        level: A nonnegative integer representing the level of SG we're 
            working with.
        address: np.array of size (level+1) representing the address 
            vector of a point in some SG graph.

    Returns:
        index: A number in the range [1, 3^(level+1)] representing the 
            TLR index of some point in level m.
    """

    # Initialize index with additional 1 added
    index = 1

    # Perform transformation similar to base 3 -> base 10
    for i in range(level+1):
        index = index + address[i] * (3 ** (level - i))

    return int(index)

def q_contract(orig_mat, base_point):
    """Computes the image of points under some contraction.

    Contracts all points represented in rows of v toward the base_point
        with a contraction ratio of 1/2.

    Args:
        orig_mat: A np.array of size (t, 2) for some integer t, with 
            each row representing the coordinate of a point
        base_point: A np.array of length 2 representing the coordinates
            of the point we're contracting toward.

    Returns:
        contract_mat: A np.array of size (t, 2) for some integer t, with 
            each row the contracted version of that in orig_mat.
    """
    # Use broadcasting to find the average points
    contract_mat = 0.5 * orig_mat + 0.5 * base_point
    return contract_mat


def sg_coordinates(level):
    """Computes the coordinates of points for SG at some level.

    Given base points q0 = [cos(5*pi/12) sin(5*pi/12)], q1=[0 0],
        q2=[cos(pi/12) sin(pi/12)], we calculate the coordinates of all
        points up to some level of the SG graph.
    
    Args:
        level: A nonnegative integer representing the level of SG we're 
            working with.

    Returns:
        coord_mat: A matrix with dimension (3^(level+1), 2), with each
            row representing coordinate of one point in V_m
    """

    # Initialize the boundary points of V_0
    q0 = np.array([np.cos(np.pi * 5 / 12), np.sin(np.pi * 5 / 12)])
    print(q0)
    q1 = np.array([0, 0])
    print(q1)
    q2 = np.array([np.cos(np.pi / 12), np.sin(np.pi / 12)])
    print(q2)

    # Stack the first three coordinates into a matrix
    coord_mat = np.stack((q0, q1, q2))

    # Each iteration produces a matrix of coordinates in the next level
    for _ in range(level):
        coord_mat = np.vstack((q_contract(coord_mat, q0),
            q_contract(coord_mat, q1), q_contract(coord_mat, q2)))

    return coord_mat


def alternate_address(level, address):
    """Computes the alternate address of a point.

    Args:
        level: A nonnegative integer representing the level of SG we're 
            working with.
        address: np.array of size (level+1) representing the address 
            vector of a point in some SG graph.
    
    Returns:
        alt_address: np.array of size (level+1) representing the 
            alternate address of the same point in some SG graph.
    
    Example:
        alternate_address(2, [0, 1, 2]) = [0, 2, 1] (F1F0q2 = F2F0q1)
        alternate_address(2, [1, 0, 0]) = [0, 1, 1] (F0F1q1 = F1F1q0)
    """

    # Make a copy of the address
    alt_address = np.copy(address)

    if (level > 0):
        # Find the index of the last address entry that is not equal 
        #   to the final value
        last_val = alt_address[level]
        temp_index = level - 1
        while (last_val == address[temp_index] and temp_index > 0):
            temp_index = temp_index - 1
        
        # If there is such an entry, interchange it with the final value
        temp = alt_address[temp_index]
        alt_address[temp_index] = alt_address[level]
        alt_address[level] = temp
    
    return alt_address 

def alternate_index(level, index):
    """Compute the alternate TLR index of a point

    Args:
        level: A nonnegative integer representing the level of SG we're 
            working with.
        index: A number in the range [1, 3^(level+1)] representing the 
            TLR index of some point in level m.
    
    Returns:
        alt_index: A number in the range [1, 3^(level+1)] representing 
            the alternate TLR index of the point in level m.
   
    Example:
        alternate_index(1, 2) = 4
        alternate_index(1, 3) = 7 
    """
    current_address = address_from_index(level, index)
    alt_address = alternate_address(level, current_address)
    alt_index = index_from_address(level, alt_address)
    return alt_index

def get_neighbors(level, address):
    """Compute the addresses of the 4 neighbors of a point

    Args:
        level: A nonnegative integer representing the level of SG we're 
            working with.
        address: np.array of size (level+1) representing the address 
            vector of a point in some SG graph.

    Returns:
        nbhd_mat: np.array of size (4, level+1) representing the 
            address of the 4 neighbors of the point. 
    """

    # Find alternate address of the point
    alt_address = alternate_address(level, address)

    # ones has all coordinates 1
    # di has only the (m+1)-th term 1, others 0
    di = np.zeros(level + 1)
    di[level] = 1

    # Initialize space for nbhd matrix
    nbhd_mat = np.zeros((4, level+1))

    # The only update happens on the last term, where q_i is replaced by 
    #   q_(i-1) and q_(i+1). Both cases reduce to mod 3 operations.
    nbhd_mat[0, :] = np.mod(address - di, 3)
    nbhd_mat[1, :] = np.mod(address + di, 3)
    nbhd_mat[2, :] = np.mod(alt_address - di, 3)
    nbhd_mat[3, :] = np.mod(alt_address + di, 3)

    return nbhd_mat

def index_alternate_level(current_level, target_level, address):
    """Finds the index of a point in another level of SG

    Args:
        current_level: A nonnegative integer representing the current
            level of SG we're working with.
        target_level: A nonnegative integer representing the target
            level of SG we're moving our point to. This should be larger
            than the current level.
        address: np.array of size (current_level+1) representing the 
            address vector of a point in some SG graph.

    Returns:
        target_index: A number in the range [1, 3^(target_level+1)] 
            representing the TLR index of the point in target_level.
    
    Example: 
        index_alternate_level(1, 1, [0, 1]) = 2
        index_alternate_level(1, 2, [0, 1]) = 5
        index_alternate_level(1, 3, [0, 1]) = 14
    """

    target_address = np.zeros(target_level+1) + address[current_level]
    target_address[0:current_level] = address[0:current_level]
    target_index = index_from_address(target_level, target_address)

    return target_index

def sg_edge_index(level, ind1, ind2):
    """Calculates an array of indices on a given edge

    Args:
        level: A nonnegative integer representing the level of SG we're 
            working with.
        ind1: a number representing the index of the first boundary 
            point (1, 2, or 3)
        ind2: a number representing the index of the second boundary
            point (1, 2, or 3, should be larger than ind1)
    
    Returns:
        index_arr: np.array of length 2^(m+1), representing array of 
            indices on the edge (each point is counted twice except 
            for the two boundary points.)
    
    Example:
        sg_edge_index(2, 0, 2) = [1, 3, 7, 9, 19, 21, 25, 27]
    """
    index_arr = np.zeros(2 ** (level+1))
    
    for j in range(2 ** (level+1)):
        temp_address = np.zeros(level+1)
        l = j
        for k in range(level+1):
            if ((l % 2) == 0):
                temp_address[level - k] = ind1
            else:
                temp_address[level - k] = ind2
            l = int(l / 2)
        index_arr[j] = index_from_address(level, temp_address)
    
    index_arr = sorted(index_arr)
    return index_arr

def rotate_address(level, address, rotate_num):
    """Rotates the address with respect to rotational symmetries of SG

    Args:
        level: A nonnegative integer representing the level of SG we're 
            working with.
        address: np.array of size (level+1) representing the address 
            vector of a point in some SG graph.
        rotate_num: A number in {0, 1, 2} representing the type of 
            rotation we're making. 0 represent no rotation, 1 represent
            counterclockwise 2*np.pi/3, 2 represent counterclockwise 
            4*np.pi/3.
    
    Returns:
        new_address: np.array of size (level+1) representing the address 
            vector of the rotated point in some SG graph.

    """
    new_address = address

    for i in range(level+1):
        new_address[i] = (address[i] + rotate_num) % 3
    return new_address



# The following methods concern the creation of polynomials
#   P_{j1} and P_{j3}. They will be moved to the monomials.py
#   file in the future.

def coeff_monomial_3_all(level, max_order):
    """Calculate the value of monomial P_{j3} up to some level
    
    Args:
        level: A nonnegative integer representing the highest level of 
            SG we would like to calculate the monomial.
        max_order: Maximal order of monomial P_{j3} that we would like 
            to calculate

    Returns:
        monomial_mat_3: np.array + size (3^{(level+1)} * (max_order + 1)
            First coordinate denote the TLR index, the second coordinate 
            denote order of polynomial.
    """
    # First initialize the case for boundary points
    q0_arr = np.zeros(max_order+1)
    q1_arr = gamma_array(max_order)
    q2_arr = -q1_arr

    # First stack the first matrix
    current_monomial_mat_3 = np.vstack((q0_arr, q1_arr, q2_arr))

    for j in range(1, level+1):
        # Make space for next-layer matrix
        temp_monomial_mat_3 = np.zeros((3 ** (j+1), max_order+1))

        # current_monomial_mat_3 (of dimension (3^j * max_order)) 
        # stores the current matrix
        for k in range(3 ** (j+1)):
            # Find the address of the point
            temp_index = k + 1
            temp_address = address_from_index(j, temp_index)
            if (k == 11):
                print(temp_address)
            # If the address is from the previous layer, we just copy
            # the value. Otherwise, we calculate the value based on 
            # the recurrent formula.
            if (temp_address[j] == temp_address[j-1]):
                # Fetch the previous index of the point
                temp_address = temp_address[0:j]
                prev_index = index_from_address(j-1, temp_address)
                
                # Set the value to be the previous value
                temp_monomial_mat_3[k, :] =  \
                    current_monomial_mat_3[prev_index-1, :]
            else:
                # Fetch the coordinates of the triangle around the point
                if (j == 1):
                    temp_address_1 = np.array([1])
                    temp_address_2 = np.array([2])
                else:
                    temp_address_1 = np.concatenate(
                        (temp_address[0:(j-1)], np.array([1])))
                    temp_address_2 = np.concatenate(
                        (temp_address[0:(j-1)], np.array([2])))

                for l in range(max_order+1):
                    # Calculate the values recursively
                    if ((temp_address[j-1] == 0 and temp_address[j] == 1) 
                            or (temp_address[j-1] == 1 and temp_address[j] == 0)):
                        # Perform the first update
                        ind_temp_address_1 = index_from_address(j-1, 
                            temp_address_1) - 1
                        temp_monomial_mat_3[k, l] = (5 ** (-l-1)) * \
                            current_monomial_mat_3[ind_temp_address_1, l]
                    elif ((temp_address[j-1] == 0 and temp_address[j] == 2) 
                            or (temp_address[j-1] == 2 and temp_address[j] == 0)):
                        # Perform the second update
                        ind_temp_address_2 = index_from_address(j-1, 
                            temp_address_2) - 1
                        temp_monomial_mat_3[k, l] = (5 ** (-l-1)) * \
                            current_monomial_mat_3[ind_temp_address_2, l]
                        if (k == 11):
                            print(temp_address_2)
                            print(ind_temp_address_2)
                            print(current_monomial_mat_3[ind_temp_address_2, l])
                    else:
                        # Perform the third update
                        rotate_add_1 = rotate_address(j-1, temp_address_2, 2)
                        rotate_add_2 = rotate_address(j-1, temp_address_1, 1)
                        ind_rotate_add_1 = index_from_address(j-1, 
                            rotate_add_1) - 1
                        ind_rotate_add_2 = index_from_address(j-1, 
                            rotate_add_2) - 1
                        temp_monomial_mat_3[k, l] = -(5 ** (-l-1)) * \
                            (current_monomial_mat_3[ind_rotate_add_1, l] + 
                            current_monomial_mat_3[ind_rotate_add_2, l])
        
        current_monomial_mat_3 = temp_monomial_mat_3

    return current_monomial_mat_3



# The following methods concern the creation of polynomials
#   P_{j1} and P_{j3}. They will be moved to the monomials.py
#   file in the future.

def coeff_monomial_1_all(level, max_order):
    """Calculate the value of monomial P_{j1} up to some level
    
    Args:
        level: A nonnegative integer representing the highest level of 
            SG we would like to calculate the monomial.
        max_order: Maximal order of monomial P_{j3} that we would like 
            to calculate

    Returns:
        monomial_mat_1: np.array + size (3^{(level+1)} * (max_order + 1)
            First coordinate denote the TLR index, the second coordinate 
            denote order of polynomial.
    """
    # TODO: Change this code to fit P_{j1}
    # First initialize the case for boundary points
    q0_arr = np.zeros(max_order+1)
    q0_arr[0] = 1
    q1_arr = alpha_array(max_order)
    q2_arr = q1_arr

    # Fetch the monomials at 3. This will be useful later.
    # stored_monomial_mat_3 = coeff_monomial_3_all(level - 1, max_order)

    # First stack the first matrix
    current_monomial_mat_1 = np.vstack((q0_arr, q1_arr, q2_arr))

    for j in range(1, level+1):
        # Make space for next-layer matrix
        temp_monomial_mat_1 = np.zeros((3 ** (j+1), max_order+1))

        # current_monomial_mat_1 (of dimension (3^j * max_order)) 
        # stores the current matrix
        for k in range(3 ** (j+1)):
            # Find the address of the point
            temp_index = k + 1
            temp_address = address_from_index(j, temp_index)

            # If the address is from the previous layer, we just copy
            # the value. Otherwise, we calculate the value based on 
            # the recurrent formula.
            if (temp_address[j] == temp_address[j-1]):
                # Fetch the previous index of the point
                temp_address = temp_address[0:j]
                prev_index = index_from_address(j-1, temp_address)
                
                # Set the value to be the previous value
                temp_monomial_mat_1[k, :] =  \
                    current_monomial_mat_1[prev_index-1, :]
            else:
                # Fetch the coordinates of the triangle around the point
                if (j == 1):
                    temp_address_1 = np.array([1])
                    temp_address_2 = np.array([2])
                else:
                    temp_address_1 = np.concatenate(
                        (temp_address[0:(j-1)], np.array([1])))
                    temp_address_2 = np.concatenate(
                        (temp_address[0:(j-1)], np.array([2])))

                for l in range(max_order+1):
                    # Calculate the values recursively
                    if (temp_address[j] == 0):
                        ind_temp_address_1 = index_from_address(j-1, 
                            temp_address_1) - 1
                        temp_monomial_mat_1[k, l] = (5 ** (-l-1)) * \
                            current_monomial_mat_1[ind_temp_address_1, l]
                    elif (temp_address[j] == 1):
                        ind_temp_address_2 = index_from_address(j-1, 
                            temp_address_2) - 1
                        temp_monomial_mat_1[k, l] = (5 ** (-l-1)) * \
                            current_monomial_mat_1[ind_temp_address_2, l]
                    else:
                        rotate_add_1 = rotate_address(j-1, temp_address_2, 2)
                        rotate_add_2 = rotate_address(j-1, temp_address_1, 1)
                        ind_rotate_add_1 = index_from_address(j-1, 
                            rotate_add_1) - 1
                        ind_rotate_add_2 = index_from_address(j-1, 
                            rotate_add_2) - 1
                        temp_monomial_mat_1[k, l] = -(5 ** (-l-1)) * \
                            (current_monomial_mat_1[ind_rotate_add_1, l] + 
                            current_monomial_mat_1[ind_rotate_add_2, l])
        
        current_monomial_mat_1 = temp_monomial_mat_1

    return current_monomial_mat_1


def progress(count, total, status=''):
    '''
    Add progress bar to function.
    
    Args:
        count: iteration of loop
        total: total number of loop iterations
        status: name of loop


    '''
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()


#print(coeff_monomial_3_all(2, 3))