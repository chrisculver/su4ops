import copy
from su4ops.constants import NS, NC
import su4ops.quark as Q
import su4ops.elemental as E
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
    'bar': True,
    'flavor': flavor,
    'color': 0,
    'spin': spin,
})


def fullVec_to_reduced(vec, mapping, NB):
  res = np.zeros(NB, dtype=np.complex64)
  for i,val in enumerate(vec):
    res[mapping[i]]+=val
  return res


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
    perm1 = perm1[:]  # copy this list so we don't mutate the original

    transCount = 0
    # Do (len - 1) transpositions
    for loc in range(len(perm0) - 1):
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
  v = copy.deepcopy(basis)
  for i in range(len(vec)):
    v[i].coef *= vec[i]
  s = ""
  for i in range(len(v)-1):
    if not np.isclose(v[i].coef, 0):
      s += "{}".format(v[i])+"+"
  if not np.isclose(v[len(v)-1].coef, 0):
    s += "{}".format(v[len(v)-1])
  print(s)

