import numpy as np
from zeta.czeta import czeta
from matrix_element import *



## basic usage

t0 = czeta()

_svec = [.1,.1,0]
_gamma = 1.1
_q2 = .9

t0.set_svec_gamma( _svec[0], _svec[1], _svec[2], _gamma)

for _l in range(5):
  for _m in range( -_l, _l+1):
    t0.set_lm( _l, _m)
    print( _l, _m, t0.evaluate(_q2))



## indexing for the M matrix

def Mindex( l, m):
  return l*(l+1) -m

def Mitolm( i):
  l = int( np.floor( np.sqrt( i)))
  m = l*(l+1)-i
  return (l, m)

lmax = 2 ## maximum angular momentum to consider
lmax2 = (lmax+1)**2 ## maximum index of matrix
M = np.zeros(( lmax2,  lmax2), dtype=complex)
for _i0 in range( lmax2):
  for _i1 in range( lmax2):
    _l0, _m0 = Mitolm( _i0)
    _l1, _m1 = Mitolm( _i1)
    M[ _i0, _i1] = matrixElem( t0, _l0, _l1, _m0, _m1, _q2)

### for consistency checks
#z00 = chop(M[0,0])
#z40 = chop(M[4,8])*7./15.
#### the following must all be the same
#print(  chop(M[4, 8])*7./15.)
#print( -chop(M[2,12])*7./(4.*np.sqrt(21.)))
#print(  chop(M[1,11])*7./(3.*np.sqrt(14.)))
#print(  chop(M[3, 9])*7./(1.*np.sqrt(210.)))

print( chop(M))



