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
    if elemental in extraBasis:
      newIdx = basis.index(extraBasis[elemental])
    else:
      newIdx = basis.index(elemental)

    newVec[newIdx] += val
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


def fullVec_to_reduced(vec, basis, extraBasis, f1=0, f2=0):
  if NC == 2:
    return su2_fullVec_to_reduced(vec, basis, extraBasis, f1, f2)
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
