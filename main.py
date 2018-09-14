import numpy as np
from zeta.czeta import czeta

print "start"
t0 = czeta()
print "test 0"

d = [0,1,1]
gam = 1.1
q2 = .9

print "test 1"
t0.set_dgam(d[0],d[1],d[2],gam)
print "test 2"

for l in range(5):
 for m in range(-l,l+1):
  t0.set_lm(l,m)
  print l,m,t0.evaluate(q2)

