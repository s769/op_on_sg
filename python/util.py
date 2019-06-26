# This file will implement helper functions written in the matlab files.
import numpy as np

# TODO

def address_from_index(level, index):
    """Computes address vector from the index.

    Returns the address vector of a point with certain TLR index in 
        a graph of SG at some fixed level.
    
    For point F_{w1} F_{w2} ... F_{wm} q_{l} (0<=wi<=2, 0<=l<=2), we can 
    represent it in two ways:
    1. Address Vector: [wm, ..., w1, l] (Note that this represents
        picking transform F_{wm}, ..., F_{w1} in this sequence and then 
        finding vertex q_{l} in the final triangle. This order is weird 
        but will prove to be useful)
    2. TLR Index: k = wm * 3^m + ... + w1 * 3 + l + 1

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
    
    For point F_{w1} F_{w2} ... F_{wm} q_{l} (0<=wi<=2, 0<=l<=2), we can 
    represent it in two ways:
    1. Address Vector: [wm, ..., w1, l] (Note that this represents
        picking transform F_{wm}, ..., F_{w1} in this sequence and then 
        finding vertex q_{l} in the final triangle)
    2. TLR Index: k = wm * 3^m + ... + w1 * 3 + l + 1

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

    return index

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
