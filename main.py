import numpy as np
from zeta.czeta import czeta
from matrix_element import *

#t0 = czeta()
#
##d = [0,1,1]
##gam = 1.1
#d = [0,0,0]
#gam = 1.0
#q2 = .9
#
#t0.set_dgam(d[0],d[1],d[2],gam)
#
#for l in range(5):
# for m in range(-l,l+1):
#  t0.set_lm(l,m)
#  print l,m,t0.evaluate(q2)

#lmax = 1
##lmax = 2
#fullString =""
#for l0 in range(0,lmax+1):
# for l2 in range(0,lmax+1):
#  for l1 in range( abs(l0-l2), l0+l2+1):
#   for m0 in range(-l0,l0+1):
#    for m1 in range(-l1,l1+1):
#     for m2 in range(-l2,l2+1):
#      print 'test',l0,l1,l2,m0,m1,m2
#      fullString = fullString +("l"+str(l0)+str(l1)+str(l2)+"m"+str(m0)+str(m1)+str(m2)+" "\
#      +("%10.6e" % wigner3j(l0,l1,l2,m0,m1,m2))+"\n")
#
#print fullString
#
#f = open('test_py_3j','w')
#f.write(fullString)
#f.close()

t0 = czeta()

d = [0,0,0]
gam = 1.0
t0.set_dgam(d[0],d[1],d[2],gam)
t0.set_lm(1,0)

q2 = .9
lmax = 2 ## maximum angular momentum to consider
lmax2 = (lmax+1)*(lmax+1) ## maximum index of matrix

## indexing for the M matrix
def Mindex(l,m):
  return l*(l+1)-m
def Mitolm(i):
  l = int( np.floor(np.sqrt(i)))
  m = l*(l+1)-i
  return (l,m)

M = np.zeros( (lmax2, lmax2), dtype=complex)
for i0 in range(lmax2):
 for i1 in range(lmax2):
  l0,m0 = Mitolm(i0)
  l1,m1 = Mitolm(i1)
  M[i0,i1] = matrixElem(t0,l0,l1,m0,m1,q2)

print M

#for l in range(5):
# for m in range(-l,l+1):
#  t0.set_lm(l,m)
#  print l,m,t0.evaluate(q2)

