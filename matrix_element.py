import cmath
import numpy as np
from scipy.special import gamma

## triangle coefficient
def triCoef(l0,l1,l2):
  return gamma(l0+l1-l2+1) *gamma(l0-l1+l2+1) *gamma(-l0+l1+l2+1) /gamma(l0+l1+l2+2)

## \sum_t (-1)^t / ( t! (l2-l1+t+m0)! (l2-l0+t-m1)! (l0+l1-l2-t)! (l0-t-m0)! (l1-t-m1)!)
## t is summed only over terms where all arguments of factorials is positive
## max t = nu can be computed before evaluation
def xSum(l0,l1,l2,m0,m1,m2):
  nu = min( l0+m0, l0-m0, l1+m1, l1-m1, l2+m2, l2-m2, l0+l1-l2, l0-l1+l2, -l0+l1+l2)
  xsum = 0.
  for t in range(nu+1): ## nu+1 terms?
    print t,(l2-l1+t+m0),(l2-l0+t-m1),(l0+l1-l2-t),(l0-t-m0),(l1-t-m1)
    xsum += np.power(-1.,t)/( \
      gamma(t+1) *gamma(l2-l1+t+m0+1) *gamma(l2-l0+t-m1+1) \
      *gamma(l0+l1-l2-t+1) *gamma(l0-t-m0+1) *gamma(l1-t-m1+1) )
  return xsum

## compute 3j from Racah formula
def wigner3j(l0,l1,l2,m0,m1,m2):
  return np.power(-1.,l0-l1-m2) \
    *np.sqrt( triCoef(l0,l1,l2) *gamma(l0+m0+1) *gamma(l0-m0+1) \
    *gamma(l1+m1+1) *gamma(l1-m1+1) *gamma(l2+m2+1) *gamma(l2-m2+1) ) \
    *xSum(l0,l1,l2,m0,m1,m2)

## clebsch-gordan coefficient calculation
def tensorC(l0,l1,l2,m0,m1,m2):
  iphase = cmath.exp( complex(0., .5*np.pi *(l0-l1+l2)) )
  return np.power(-1.,m2) *np.sqrt( (2.*l0+1) *(2.*l1+1) *(2.*l2+1) ) \
    *wigner3j(l0,l1,l2,0,0,0) *wigner3j(l0,l1,l2,m0,m1,m2) *iphase

## computation of matrix elements
## cz is a czeta class implementation
def matrixElem(cz,l0,l2,m0,m2,q2):
  pref = np.power(-1.,l0) /np.power(np.pi,1.5)
  msum = complex(0.,0.)
  for l1 in range( abs(l0-l2), l0+l2+1):
    for m1 in range( -l1, l1+1):
      cz.set_lm(l1,m1)
      zeta = cz.evaluate(q2)
      tc = tensorC(l0,l1,l2,m0,m1,m2)
      iphase = cmath.exp( complex(0., .5*np.pi*l1) )
      msum += np.power( q2, -.5*(l1+1.)) *zeta *tc *pref *iphase
  return msum

