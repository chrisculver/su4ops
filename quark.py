import math
import numpy as np
from constants import NS, gammas, NF, flavorLabels, colorLabels


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
      if len(self.quarks) == 3:
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


class Quark():
    def __init__(self, p):
        if not isinstance(p, dict):
            raise TypeError('Need to pass a dictionary of quark properties')

        self.barred = p['bar']
        self.flavor = p['flavor']
        self.color = p['color']
        self.spin = p['spin']

    def spatial_rotate(self, gElem):
        vec = [1 if i == self.spin else 0 for i in range(NS)]
        if self.barred:
            rotVec = np.matmul(np.transpose(vec), Gamma(gElem, self.barred))
        else:
            rotVec = np.matmul(Gamma(gElem, self.barred), vec)

        return rotVec

    def g_parity_rotate(self):
      spinVec = [1 if i == self.spin else 0 for i in range(NS)]
      spinVec = np.matmul(gammas[2], spinVec)

      flavorVec = []
      if NF == 2:
        flavorVec = [0 if f == self.flavor else -1 for f in [0, 1]]
      else:
        raise ValueError("Gparity not implemented for NF={}".format(NF))

      antiVec = [0 if i == self.barred else 1 for i in [
          True, False]]  # swap quark w/ anti-quark

      vec = np.kron(spinVec, flavorVec)
      vec = np.kron(vec, antiVec)

      return vec

    def __eq__(self, other):
      return self.barred == other.barred and self.flavor == other.flavor and self.color == other.color and self.spin == other.spin

    def __hash__(self):
      return hash((self.barred, self.flavor, self.color, self.spin))

    def __str__(self):
        s = ''
        if self.barred:
            s += '\\bar{'
        s += flavorLabels[self.flavor]
        if self.barred:
            s += '}'
        s += '_'+'{' + colorLabels[self.color] + ',' + str(self.spin) + '}'

        return s


def Gamma(gElem, barred=False):
    if gElem.identifier['parity'] in [1, None]:
        return GammaProper(gElem, barred)
    else:
        return np.matmul(gammas[4], GammaProper(gElem, barred))


def GammaProper(gElem, barred=False):
    id = np.identity(NS)
    t0 = id*math.cos(gElem.identifier['angle']/2.)

    dir = gElem.identifier['direction']
    dir = dir/np.linalg.norm(dir)
    gammaAngle = dir[0]*np.matmul(gammas[2], gammas[3])+dir[1]*np.matmul(
        gammas[3], gammas[1])+dir[2]*np.matmul(gammas[1], gammas[2])

    t1 = gammaAngle*math.sin(gElem.identifier['angle']/2.)

    if barred:
        return t0-t1

    return t0+t1
