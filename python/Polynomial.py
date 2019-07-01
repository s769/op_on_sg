import numpy as np
from innerprods import lis2str, inner_dict, symmetrize, vals_dict, norm_dict
import sympy as sp
from sympy import Rational as Rat
import tqdm
'''
This file contains the class Polynomial which is used to create orthogonal 
    polynomials with respect to various inner products.
'''


class Polynomial:
    # This stores the Gram Matrices for various inner products

    GM = {}

    @staticmethod
    def has_GM(n, lam=np.array([1])):
        '''
        This function checks whether a Gram Matrix of a given size has 
            been built for a given inner product.
        
        Args:
            n: size of Gram Matrix required
            lam: array of lambda values for the inner product (default 
                is [1] which is the regular Sobolev inner product).
        Returns:
            Boolean value expressing whether such a Gram Matrix has been 
                constructed.
        '''
        try:
            return (not Polynomial.GM is None) and\
            Polynomial.GM[lis2str(lam)].shape[0] >= n
        except KeyError:
            return False

    def __init__(self, coefs, j, k, lam=np.array([1])):
        '''
        Constructor for Polynomial.

        Args:
            coefs: np.array of coefficients of the polynomial in terms of 
                the monomial basis P_jk. (Note: the basis is ordered 
                {P_01, P_02, P_03, P_11, P_12, P_13, ...}).
            j: degree of polynomial. If the length of coefs does not 
                equal 3*j + 3, zeros are added to the end of coefs.
            k: family of polynomial. Since most of the polynomials we 
                will deal with will come only from a certain family 
                (k = 1, 2, or 3), setting the k value will allow us to 
                deal with only the basis {P_0k, P_1k, ...}.
            lam: np.array of lambda_values used for calculating inner 
                products. If lam contains only one value, it is set as 
                a scalar. The default value is 1 (corresponding to the 
                regular Sobolev inner product).
        '''
        if not len(coefs) == 3*j+3:
            ad = sp.zeros(3*j+3 - len(coefs), 1)
            coefs.col_join(ad)
        self.coefs = coefs
        self.j = j
        self.k = k
        if len(lam) == 1:
            lam = lam[0]
        self.lam = lam

    def __add__(self, other):
        '''
        Adds 2 polynomials by adding their coefficients
        '''
        return Polynomial(self.coefs+other.coefs, self.j, self.k, self.lam)

    @staticmethod
    def basis_inner(j, i, k, ip, lam=np.array([1])):
        '''
        Computes any Sobolev inner product of <P_ji, P_ki'>.
        Args:
            j, i, k, ip: correspond to j, i, k, i' in P_ji, P_ki' 
                respectively.
            lam: np.array of lambda values for the generalized Sobolev 
                inner product. The default value is 1 (corresponding to 
                the regular Sobolev inner product). 
                If lam = np.array([0]), this is the L2 inner product.

        Returns:
            Sobolev inner product of <P_ji, P_ki'> with given lambda 
                values.
        '''

        # inner_dict is found in innerprods.py
        try:
            inner_func = inner_dict[(i, ip)] #Example: if i=1 and ip = 2 then inner_func is <P_{j1}, P_{k2}>_{L2} 
        except KeyError:
            inner_func = inner_dict[(ip, i)]

        res = inner_func(j, k) #Initialize to standard L2 inner product
        for i in range(1, len(lam)+1): 
            res += inner_func(j-i, k-i)*lam[i-1] #This is the step where you add the sobolev terms

        return res

    @staticmethod
    def slow_inner(arr1, arr2, lam=np.array([1])):
        '''
        Computes any inner product between two polynomials represented 
            as arrays of coefficients.

        Args:
            arr1, arr2: arrays of coefficients of polynomials in monomial 
                basis.(Note: the basis is ordered 
                {P_01, P_02, P_03, P_11, P_12, P_13, ...}).
            lam: np.array of lambda values for the generalized Sobolev 
                inner product. The default value is 1 (corresponding to 
                the regular Sobolev inner product). 
                If lam = np.array([0]), this is the L2 inner product.

        Returns:
            Sobolev inner product between the polynomials represented 
                by arr1, and arr2 with the given lambda values.
        '''
        res = 0
        for ind1 in range(len(arr1)):
            for ind2 in range(len(arr2)):
                if arr1[ind1] == 0 or arr2[ind2] == 0:
                    continue
                #Here j,k,i, and ip are recovered based on the way the array of co-efficients indexes the polynomials
                j = int(np.floor(ind1/3)) 
                k = int(np.floor(ind2/3))
                i = int(ind1 % 3 + 1)
                ip = int(ind2 % 3 + 1)
                res += \
                    arr1[ind1]*arr2[ind2]\
                    * Polynomial.basis_inner(j, i, k, ip, lam)

        return res

    @staticmethod
    def build_GM(n, lam=np.array([1])):
        '''
        Constructs Gram Matrix of a given size for a given generalized 
            Sobolev inner product. This is used to later compute inner 
            products more quickly.
        The Gram Matrix is stored in the dictionary Polynomial.GM keyed 
            by the string representing the values of lambda.

        Args:
            n: size of Gram Matrix
            lam: np.array of lambda values for the generalized Sobolev 
                inner product. The default value is 1 (corresponding to 
                the regular Sobolev inner product).
                If lam = np.array([0]), this is the L2 inner product.
        '''
        if Polynomial.has_GM(3*n+3, lam):
            return


        arr = sp.zeros(3*n+3, 3*n+3)
        #The following if else statements make an array lam_arr that represents the weights of all the integrals 
        # in the inner product formula. This is required because lam only describes the weights on the integrals with positive order laplacians because
        # integrals with positive order laplacians as we consider the weight on the L2 inner product to be 1. 

        if not (np.array_equal(lam, np.array([1])) \
                or np.array_equal(lam, np.array([0]))):
            lam_arr = [np.array([1]), np.array([0]), lam]
        else:
            lam_arr = [np.array([1]), np.array([0])]

        for lamb in lam_arr:
            for ind1 in range(n):
                for ind2 in range(n):
                    if ind1 <= ind2:
                        j = int(np.floor(ind1/3))
                        k = int(np.floor(ind2/3))
                        i = int(ind1 % 3 + 1)
                        ip = int(ind2 % 3 + 1)
                        arr[ind1, ind2] = \
                            Polynomial.basis_inner(j, i, k, ip, lamb)
            Polynomial.GM[lis2str(lamb)] = symmetrize(arr)

      
        return

    @staticmethod
    def build_condensed_GM(n, i, lam=np.array([1])):
        '''
        When we work with only polynomials from a certain family 
            (k = 1, 2, or 3), it is more convenient to only work with 
            the basis {P_0k, P_1k,...}. This function creates the 
            Gram Matrix for a given generalized Sobolev inner product.

        Args:
            n: size of Gram Matrix required
            i: family of polynomials (i represents k in the preceding 
                paragraph)
            lam: np.array of lambda values for the generalized Sobolev 
                inner product. The default value is 1 (corresponding to 
                the regular Sobolev inner product).
                If lam = np.array([0]), this is the L2 inner product.
        '''
        if Polynomial.has_GM(3*n+3, lam):
            return

        arr = sp.zeros(n, n)
        if not (np.array_equal(lam, np.array([1])) \
                or np.array_equal(lam, np.array([0]))):
            lam_arr = [np.array([1]), np.array([0]), lam]
        else:
            lam_arr = [np.array([1]), np.array([0])]
        for lamb in lam_arr:
            for ind1 in tqdm.tqdm(range(n)):
                for ind2 in range(n):
                    if ind1 <= ind2:
                        arr[ind1, ind2] = \
                            Polynomial.basis_inner(ind1, i, ind2, i, lamb)

            Polynomial.GM[lis2str(lamb)] = symmetrize(arr)


    @staticmethod
    def fast_inner(arr1, arr2, GM):
        '''
        Computes any inner product between two polynomials represented 
            as arrays of coefficients.

        Args:
            arr1, arr2: arrays of coefficients of polynomials in monomial 
                basis.
            GM: Gram Matrix of inner product (Note: the baisis order 
                for arr1 and arr2 must match that of GM)

        Returns:
            Sobolev inner product between the polynomials represented 
                by arr1, and arr2 for the inner product represented by GM.
        '''
        return arr1.T @ (GM @ arr2)

    # This function pads the coefs arrays of 2 objects with zeros so that they have the same length
    @staticmethod
    def pad(obj1, obj2):
        if len(obj1.coefs) > len(obj2.coefs):
            obj2.coefs.col_join(sp.zeros(len(obj1.coefs)-len(obj2.coefs)))
        elif len(obj1.coefs) < len(obj2.coefs):
            obj1.coefs.col_join(sp.zeros(len(obj2.coefs)-len(obj1.coefs)))
        arr1 = obj1.coefs
        arr2 = obj2.coefs

        return arr1, arr2

    def inner(self, other):
        '''
        Computes the inner product between 2 Polynomial objects.
            If the Gram Matrix is available, this function uses 
            fast_inner. Otherwise, this function uses slow_inner.    

        '''
        arr1, arr2 = Polynomial.pad(self, other)

        n = len(arr1)
        GM_status = Polynomial.has_GM(n, self.lam)
        if GM_status:
            print('Using Fast Inner Product')
            GM = Polynomial.GM[lis2str(self.lam)][:n, :n]
            return Polynomial.fast_inner(arr1, arr2, GM)
        print('Using Slow Inner Product')
        return Polynomial.slow_inner(arr1, arr2, self.lam)

    # This function computes the norm of a Polynomial object.
    def norm(self):
        return sp.sqrt(self.inner(self))

    def get_condensed_coefs(self):
        '''
        This function returns the array of coefficients with respect to 
            the condensed basis {P_0k, P_1k,... } when we work with 
            polynomials from only one family (k = 1, 2, or 3).

        '''
        return self.coefs[self.k-1::3]

    def value(self, j, i=1):
        '''
        This function calculates the value of a polynomial on a 
            boundary point of SG

        Args:
            j: number of polynomials to evaluate
            i: index of boundary point i = 0, 1, 2

        Returns:
            value of the polynomial at q_i

        '''
        # array, vals_dict are found in innerprods.py
        ccoefs = self.get_condensed_coefs()
        valarr = np.array([vals_dict[self.k](m) for m in range(j+1)])
        res = ccoefs.dot(valarr)
        if i == 0:
            return ccoefs[0] if self.k == 0 else 0
        if i == 1:
            return res
        if i == 2:
            return res if self.k == 0 or self.k == 1 else -res

    def dnvalue(self, j, i=1):
        '''
        This function calculates the value of the normal derivative of 
            a polynomial on a boundary point of SG

        Args:
            j: number of polynomials to evaluate
            i: index of boundary point i = 0, 1, 2

        Returns:
            value of the normal derivative of polynomial at q_i

        '''
        # array, norm_dict are found in innerprods.py
        ccoefs = self.get_condensed_coefs()
        dnarr = np.array([norm_dict[self.k](m) for m in range(j+1)])
        res = ccoefs.dot(dnarr)
        if i == 0:
            return ccoefs[0] if self.k == 1 else 0
        if i == 1:
            return res
        if i == 2:
            return res if self.k == 0 or self.k == 1 else -res

    @staticmethod
    def get_condensed_GM(k, lam=np.array([1])):
        '''
        This function takes the existing Gram Matrix for a given 
            generalized Sobolev inner product expressed with respect to 
            the orginal monomial basis 
            {P_01, P_02, P_03, P_11, P_12, P_13, ...} and 
            returns the condensed matrix expressed with respect 
            to the basis {P_0k, P_1k,... }.

        Args:
            k: index corresponding to 'k' in the preceding paragraph

            lam: np.array of lambda values for the generalized Sobolev 
                inner product. The default value is 1 (corresponding to 
                the regular Sobolev inner product).
                If lam = np.array([0]), this is the L2 inner product.

        Returns:
            Gram Matrix with of given inner product with respect 
                to the basis {P_0k, P_1k,... }.
        '''
        GM = Polynomial.GM[lis2str(lam)]
        return GM[k-1::3, k-1::3]
