import cmath
import numpy as np
from scipy.special import gamma
from sympy.physics.wigner import gaunt



## if abs2 is less than some tolerance, return 0. of the same type
def isSmall(x):
  if np.abs(x) < 1e-8:
    return True
  return False

## if abs2 is less than some tolerance, return 0. of the same type
def chop(x):
  if isinstance(x,np.ndarray):
    return np.array([ chop(y) for y in x ])
  if isinstance(x,complex):
    rx = ( 0. if isSmall(np.real(x)) else np.real(x) )
    ix = ( 0. if isSmall(np.imag(x)) else np.imag(x) )
    return complex(rx,ix)
  else:
    if isSmall(x):
      return x*0.
  return x

## clebsch-gordan coefficient calculation
## gaunt is \propto WigJ * WigJ
def tensorC(l0,l1,l2,m0,m1,m2):
  iphase = np.sqrt(4.*np.pi)*np.power(-1.,m2)* cmath.exp( complex(0., .5*np.pi *(l0-l1+l2)) )
  return gaunt(l0,l1,l2,m0,m1,m2) *iphase

## computation of matrix elements
## cz is a czeta class implementation
def matrixElem(cz,l0,l2,m0,m2,q2):
  pref = np.power(-1.,l0) /np.power(np.pi,1.5)
  msum = complex(0.,0.)
  for l1 in range( abs(l0-l2), l0+l2+1):
    for m1 in range( -l1, l1+1):
      cz.set_lm(l1,m1)
      zeta = chop(cz.evaluate(q2))
      if isSmall(zeta): ## don't deal with terms that are unphysical
        continue
      tc = tensorC(l0,l1,l2,m0,m1,-m2) ## NB: -1.*m2
      iphase = cmath.exp( complex(0., .5*np.pi*l1) )
      msum += np.power( q2, -.5*(l1+1.)) *zeta *tc *pref *iphase
  return msum

