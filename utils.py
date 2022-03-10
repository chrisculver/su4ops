import copy
from constants import NS, NC
import quark as Q
import numpy as np
import sympy as sp
# from https://www.bernardosulzbach.com/heaps-algorithm/


def _swap(elements, i, j):
    elements[i], elements[j] = elements[j], elements[i]

# from https://www.bernardosulzbach.com/heaps-algorithm/


def _generate_permutations(elements, n):
    # As by Robert Sedgewick in Permutation Generation Methods
    c = [0] * n
    yield elements
    i = 0
    while i < n:
        if c[i] < i:
            if i % 2 == 0:
                _swap(elements, 0, i)
            else:
                _swap(elements, c[i], i)
            yield elements
            c[i] += 1
            i = 0
        else:
            c[i] = 0
            i += 1


def permutations(elems):
    elements = copy.deepcopy(elems)  # leave original list in order
    return _generate_permutations(elements, len(elements))


def quark(spin, flavor=0): return Q.Quark({
    'bar': False,
    'flavor': flavor,
    'color': 0,
    'spin': spin,
})


def su2_fullVec_to_reduced(vec, basis, extraBasis, f1=0, f2=0):
  newVec = [0 for b in basis]
  for i, val in enumerate(vec):
    s0 = i % NS
    tmp = i//NS
    s1 = tmp % NS
    elemental = Q.Elemental(1, [quark(s0, f1), quark(s1, f2)])
    newIdx = 0
    sign = 1
    if elemental in extraBasis:
      newIdx = basis.index(extraBasis[elemental]['elemental'])
      sign = extraBasis[elemental]['sign']
    else:
      newIdx = basis.index(elemental)

    newVec[newIdx] += sign*val
  return np.array(newVec).round(8)

def su3_fullVec_to_reduced(vec, basis, extraBasis, f1=0, f2=0, f3=0):
  newVec = [0 for b in basis]
  for i, val in enumerate(vec):
    s0 = i % NS
    tmp = i//NS
    s1 = tmp % NS
    tmp = tmp//NS
    s2 = tmp % NS

    elemental = Q.Elemental(1, [quark(s0,f1), quark(s1,f2), quark(s2,f3)])
    newIdx = 0
    sign = 1
    if elemental in extraBasis:
      newIdx = basis.index(extraBasis[elemental]['elemental'])
      sign = extraBasis[elemental]['sign']
    else:
      newIdx = basis.index(elemental)

    newVec[newIdx] += sign*val
  return np.array(newVec).round(8)


def su4_fullVec_to_reduced(vec, basis, extraBasis):
  newVec = [0 for b in basis]
  for i, val in enumerate(vec):
    s0 = i % NS
    tmp = i//NS
    s1 = tmp % NS
    tmp = tmp//NS
    s2 = tmp % NS
    tmp = tmp//NS
    s3 = tmp % NS

    elemental = Q.Elemental(1, [quark(s0), quark(s1), quark(s2), quark(s3)])
    newIdx = 0
    if elemental in extraBasis:
      newIdx = basis.index(extraBasis[elemental])
    else:
      newIdx = basis.index(elemental)

    #TODO: I think this needs to be multiplied by +/- 1 depending on Grassman?
    newVec[newIdx] += val
  return np.array(newVec).round(8)


<<<<<<< HEAD
def fullVec_to_reduced(vec, basis, extraBasis,f1=0,f2=0,f3=0):
  if NC == 2:
    return su2_fullVec_to_reduced(vec, basis, extraBasis,f1,f2)
  elif NC == 3:
    return su3_fullVec_to_reduced(vec, basis, extraBasis,f1,f2,f3)
=======
def fullVec_to_reduced(vec, basis, extraBasis, f1=0, f2=0):
  if NC == 2:
    return su2_fullVec_to_reduced(vec, basis, extraBasis, f1, f2)
>>>>>>> c054202b26b4147517d88905e2d4865e2d0186f6
  elif NC == 4:
    return su4_fullVec_to_reduced(vec, basis, extraBasis)
  else:
    raise ValueError("Reducing NC={} vector not implemented".format(NC))
# [A,B],
# A -> A + B - C,  [1,1,-1]
# C === B
# [1,0]


def makeRepMat(basis, extraBasis, gElem, id, f1=0, f2=0):
  rotMat = []
  refMat = []
  for b in basis:
    rotMat.append(fullVec_to_reduced(
      b.spatial_rotate(gElem), basis, extraBasis, f1, f2))
    refMat.append(fullVec_to_reduced(
        b.spatial_rotate(id), basis, extraBasis, f1, f2))

  res = np.zeros((len(basis), len(basis)), dtype=complex)
  for r in range(len(basis)):
    for c in range(len(basis)):
      res[r, c] = np.dot(rotMat[r], refMat[c])

  return np.transpose(res)


def projectorMat(irrep, rep, group, row=0):
    dim = len(rep[0])
    res = np.zeros((dim, dim), dtype=complex)
    for i, elem in enumerate(group.elements):
        #print(type(float(elem.irreps[irrep][row,row])))
        res += complex(elem.irreps[irrep][row, row])*np.transpose(rep[i])
    return (res*len(irrep)/len(group.elements)).round(8)


def operators(irrep, rep, group, row=0):
  proj = projectorMat(irrep, rep, group, row)
  rrefMat = sp.Matrix(proj).rref()[0].tolist()
  for row in rrefMat[:]:
    if (row == [0 for i in range(len(row))]):
      rrefMat.remove(row)
  rrefMat = np.array(rrefMat).astype(complex)
  return rrefMat


def op_basis_map(op):
  basis_map = {}
  for i, val in enumerate(op):
    if not np.isclose(val, 0):
      basis_map[i] = val
  return basis_map




# from https://stackoverflow.com/questions/1503072/how-to-check-if-permutations-have-equal-parity #you'll never guess what my google search was to find this.
def arePermsEqualParity(perm0, perm1):
    """Check if 2 permutations are of equal parity.

    Assume that both permutation lists are of equal length
    and have the same elements. No need to check for these
    conditions.

    :param perm0: A list.
    :param perm1: Another list with same elements.

    :return: True if even parity, False if odd parity.
    """
    perm1 = perm1[:] ## copy this list so we don't mutate the original

    transCount = 0
    for loc in range(len(perm0) - 1):                         # Do (len - 1) transpositions
        p0 = perm0[loc]
        p1 = perm1[loc]
        if p0 != p1:
            sloc = perm1[loc:].index(p0)+loc          # Find position in perm1
            perm1[loc], perm1[sloc] = p0, p1          # Swap in perm1
            transCount += 1

    # Even number of transpositions means equal parity
    if (transCount % 2) == 0:
        return True
    else:
        return False

def print_vec(vec, basis):
  s=""
  for i in range(len(vec)-1):
    if not np.isclose(vec[i],0):
      s+="{}*{}".format(vec[i],basis[i])+"+"
  s+="{}*{}".format(vec[i],basis[i])
  print(s)
