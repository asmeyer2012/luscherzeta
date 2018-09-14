import numpy as np
from zeta.czeta import czeta
from matrix_element import *

t0 = czeta()

d = [0,1,1]
gam = 1.1
q2 = .9

t0.set_dgam(d[0],d[1],d[2],gam)

for l in range(5):
 for m in range(-l,l+1):
  t0.set_lm(l,m)
  print l,m,t0.evaluate(q2)

