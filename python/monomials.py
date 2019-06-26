from recursions import *



@mem
def big_recursion(j):
    if j == 0: return np.array([7/45, 4/45, 2/5, 1/5])
    
    res3 = 0
    res4 = 0
    vec2 = np.zeros(2)
    for l in range(j):
        p, q = big_recursion(l)[-2:]
        a, b = big_recursion(j-1-l)[:2]
        res3 += (4*a + 3*b)*p + (a + 2*b)*q
        res4 += (2*a + 4*b)*p + (3*a + b)*q
    b = big_recursion(j-1)[1]
    vec2[0] = -2/5*b - 1/5*res3
    vec2[1] = -1/5*b - 1/5*res4



    vec = np.zeros(2)
    res1 = 0
    res2 = 0

    for l in range(j):
        if l == 0: 
            p, q = vec2
        else: 
            p, q = big_recursion(j-l)[-2:]
        a, b  = big_recursion(l)[:2]
        res1 += (2*p + q)*(a+2*b)
        res2 += (p + 2*q)*(a + 2*b)
    vec[0] = 2/(3*(5**j - 1))*res1
    vec[1] = 10/(3*(5**(j+1) - 1))*res2

    coefs = np.array([[1, 2], [-4, 7]])
    ab_arr = np.linalg.inv(coefs)@vec

    

    return np.append(ab_arr, vec2)
# for j in range(7):
#     print(big_recursion(j))