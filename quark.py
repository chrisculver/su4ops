import scipy
import gammas
import math
import numpy as np
from constants import *


class Expression():
    def __init__(self, elementals):
        self.elementals = elementals

    def __str__(self):
        s = ''
        for e in self.elementals[:-1]:
            s+='{}+'.format(e)
        s+='{}'.format(self.elementals[-1])
        return s

    def __add__(self, other):
        if isinstance(other,int) or isinstance(other,float):
            if math.isclose(other,0):
                return self
            else:
                raise ValueError('Trying to add {} to expression'.format(other))

        if isinstance(other, Elemental):
            if other in self.elementals:
                idx = self.elementals.index(other)
                self.elementals[idx].coef+=other.coef
                if math.isclose(self.elementals[idx].coef,0):
                    self.elementals.pop(idx)
                return self

        if isinstance(other, Expression):
            for e in other.elementals:
                self=self+e
            return self




class Elemental():
    def __init__(self, coef, quarks):
        self.coef = coef
        self.quarks = quarks

    def spatial_rotate(self, gElem):
        rotationVecs = []
        for q in self.quarks:
            rotationVecs.append(q.spatial_rotate(gElem))

        res = rotationVecs[0]
        for i in range(1,len(rotationVecs)):
            res = np.kron(res,rotationVecs[i])

        #res = rotationVecs[len(rotationVecs)-1]
        #for i in range(len(rotationVecs)-2,-1,-1):
        #    res=np.kron(rotationVecs[i],res)

        return res

    def __str__(self):
        s = '{}*'.format(self.coef)
        for q in self.quarks[:-1]:
            s+='{}*'.format(q)
        s+='{}'.format(self.quarks[-1])
        return s

    def __rmul__(self, other):
        if isinstance(other, float):
            if math.isclose(self.coef,0):
                return 0
            else:
                return self

    def __add__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            if math.isclose(other,0):
                return self
            else:
                raise ValueError('Adding constant {} to elemental'.format(other,self))

        if isinstance(other, Elemental):
            return Expression([self,other])


    __radd__ = __add__





class Quark():
    def __init__(self, p):
        if not isinstance(p, dict):
            raise TypeError('Need to pass a dictionary of quark properties')

        self.barred = p['bar']
        self.flavor = p['flavor']
        self.color = p['color']
        self.spin = p['spin']


    def spatial_rotate(self, gElem):
        vec = [1 if i==self.spin else 0 for i in range(NS)]
        if self.barred:
            rotVec = np.matmul(np.transpose(vec),Gamma(gElem,self.barred))
        else:
            rotVec = np.matmul(Gamma(gElem,self.barred), vec)


        return rotVec


    def __str__(self):
        s=''
        if self.barred:
            s+='\\bar{'
        s+=flavorLabels[self.flavor]
        if self.barred:
            s+='}'
        s+='_'+'{' + colorLabels[self.color] + ',' + str(self.spin) + '}'

        return s





def Gamma(gElem, barred=False):
    if gElem.identifier['parity'] in [1,None]:
        return GammaProper(gElem,barred)
    else:
        return np.matmul(gammas[4],GammaProper(gElem,barred))


def GammaProper(gElem, barred=False):
    id = np.identity(NS)
    t0 = id*math.cos(gElem.identifier['angle']/2.)

    dir = gElem.identifier['direction']
    dir = dir/np.linalg.norm(dir)
    gammaAngle = dir[0]*np.matmul(gammas[2],gammas[3])+dir[1]*np.matmul(gammas[3],gammas[1])+dir[2]*np.matmul(gammas[1],gammas[2])

    t1 = gammaAngle*math.sin(gElem.identifier['angle']/2.)

    if barred:
        return t0-t1

    return t0+t1
