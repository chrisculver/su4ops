import numpy as np
from su4ops.constants import NS


class Elemental():
    def __init__(self, coef, quarks):
        self.coef = coef
        self.quarks = quarks

    def spatial_rotate(self, gElem):
        rotationVecs = []
        for q in self.quarks:
            rotationVecs.append(q.spatial_rotate(gElem))

        res = rotationVecs[0]
        for i in range(1, len(rotationVecs)):
            #this choice of ordering for the kronecker products means, for mesons
            # the basis elements are phi_00, phi_10, phi_20 ..., aka row-col ordering
            res = np.kron(rotationVecs[i], res)

        return res

    def g_parity_rotate(self):
      rotationVecs = []
      for q in self.quarks:
          rotationVecs.append(q.g_parity_rotate())

      res = rotationVecs[0]
      for i in range(1, len(rotationVecs)):
          #this choice of ordering for the kronecker products means, for mesons
          # the basis elements are phi_00, phi_10, phi_20 ..., aka row-col ordering
          res = np.kron(rotationVecs[i], res)

      return res

    def get_index(self):
      idx = 0
      for i in range(len(self.quarks)-1, -1, -1):
        idx += (NS**i)*self.quarks[i].spin
      return idx

    def __hash__(self):
      if len(self.quarks) == 2:
        return hash((self.coef, self.quarks[0], self.quarks[1]))
      if len(self.quarks) == 3:
        return hash((self.coef, self.quarks[0], self.quarks[1], self.quarks[2]))
      if len(self.quarks) == 4:
        return hash((self.coef, self.quarks[0], self.quarks[1], self.quarks[2], self.quarks[3]))

    def __eq__(self, other):
      return self.coef == other.coef and self.quarks == other.quarks

    def __str__(self):
        s = '{}*'.format(self.coef)
        for q in self.quarks[:-1]:
            s += '{}*'.format(q)
        s += '{}'.format(self.quarks[-1])
        return s
