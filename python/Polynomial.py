from innerprods import *

class Polynomial:
  
  GM = None


    
  def has_GM(n):
    return (not Polynomial.GM is None) and\
        Polynomial.GM['sob'].shape[0] >= n
  
  def __init__(self, coefs, j, k, lam=np.array([1])):
    if not len(coefs) == 3*j+3: 
      ad = np.zeros(3*j+3 - len(coefs))
      coefs = np.append(coefs, ad)
    self.coefs = coefs
    self.j = j
    self.k = k
    if len(lam) == 1: lam = lam[0]     
    self.lam = lam
    
  def __add__(self, other):
    return Polynomial(self.coefs+other.coefs, self.j)
  
  


  def basis_inner(j, i, k, ip, lam=np.array([1])):
    try:
      inner_func = inner_dict[(i, ip)]
    except KeyError:
      inner_func = inner_dict[(ip, i)]
      
    res = inner_func(j, k)
    for i in range(1, len(lam)+1):
      res += inner_func(j-i, k-i)*lam[i-1]
      
    return res
  
  def slow_inner(arr1, arr2, lam=np.array([1])):
    res = 0
    for ind1 in range(len(arr1)):
      for ind2 in range(len(arr2)):
        if arr1[ind1] == 0 or arr2[ind2] == 0: continue
        j = int(np.floor(ind1/3))
        k = int(np.floor(ind2/3))
        i = int(ind1%3 + 1)
        ip = int(ind2%3 + 1)
        res += \
        arr1[ind1]*arr2[ind2]\
        *Polynomial.basis_inner(j, i, k, ip, inner)
        
    return res
  
  def build_GM(n):
    if Polynomial.has_GM(n): return
    GM = {}
    
    arr = np.zeros((3*n+3, 3*n+3))
    
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
            i = int(ind1%3 + 1)
            ip = int(ind2%3 + 1)
            arr[ind1, ind2] = \
              Polynomial.basis_inner(j, i, k, ip, lamb)
      GM[lis2int(lamb)] = symmetrize(arr)
    
    Polynomial.GM = GM
    return
  
  def build_condensed_GM(n, i, lam=np.array([1])):
    GM = {}
    arr = np.zeros((n, n))
    if not (np.array_equal(lam, np.array([1])) \
            or np.array_equal(lam, np.array([0]))):
      lam_arr = [np.array([1]), np.array([0]), lam]
    else:
      lam_arr = [np.array([1]), np.array([0])]
    for lamb in lam_arr:
      for ind1 in range(n):
        for ind2 in range(n):
          if ind1 <= ind2:
            arr[ind1, ind2] = \
              Polynomial.basis_inner(ind1, i, ind2, i, lamb)
            
      GM[lis2int(lamb)] = symmetrize(arr)
      
    Polynomial.GM = GM
  
  def fast_inner(arr1, arr2, GM):
    return arr1.T @ GM @ arr2
  
  
  def pad(obj1, obj2):
    if len(obj1.coefs) > len(obj2.coefs):
      arr1 = obj1.coefs
      arr2 = np.append(obj2.coefs,\
                       np.zeros(len(obj1.coefs)-len(obj2.coefs)))
    elif len(obj1.coefs) < len(obj2.coefs):
      arr1 = np.append(obj1.coefs,\
                       np.zeros(len(obj2.coefs)-len(obj1.coefs)))
      arr2 = obj2.coefs
    else:
      arr1 = obj1.coefs
      arr2 = obj2.coefs
      
    return arr1, arr2
  
  def inner(self, other):
    arr1, arr2 = Polynomial.pad(self, other)
    
    n = len(arr1)
    GM_status = Polynomial.has_GM(n)
    print(GM_status)
    if GM_status:
      print('Using Fast Inner Product')
      GM = Polynomial.GM[lis2int(self.lam)][:n, :n]
      return Polynomial.fast_inner(arr1, arr2, GM)
    print('Using Slow Inner Product')
    return Polynomial.slow_inner(arr1, arr2, self.lam)
  
  def norm(self):
    return np.sqrt(inner(self, self))
  
  
  def get_condensed_coefs(self):
    return self.coefs[self.k-1::3]
  
  def value(self, i=1):
    ccoefs = self.get_condensed_coefs()
    valarr = array([vals_dict[self.k](m) for m in range(j+1)])
    res = ccoefs.dot(valarr)
    if i == 0:
      return np.sum(ccoefs) if self.k == 0 else 0
    if i == 1:
      return res
    if i == 2:
      return res if self.k == 0 or self.k == 1 else -res
  
  def dnvalue(self, i=1):
    ccoefs = self.get_condensed_coefs()
    dnarr = array([norm_dict[self.k](m) for m in range(j+1)])
    res = ccoefs.dot(dnarr)
    if i == 0:
      return np.sum(ccoefs) if self.k == 1 else 0
    if i == 1:
      return res
    if i == 2:
      return res if self.k == 0 or self.k == 1 else -res
  
  def get_condensed_GM(k, lam=np.array([1])):
    GM = Polynomial.GM[lis2int(lam)]
    return GM[k-1::3, k-1::3]
      