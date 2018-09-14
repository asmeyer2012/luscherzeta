import ctypes
from ctypes import c_int,c_double,c_void_p
#import numpy as np

lib = ctypes.cdll.LoadLibrary('./zeta/libczeta.so')

lib.czeta.argtypes = [ ]
lib.set_dgam.argtypes = [ c_void_p, c_int, c_int, c_int, c_double]
lib.set_lm.argtypes = [ c_void_p, c_int, c_int]
lib.zeta_eval.argtypes = [ c_void_p, c_double]
## uppercase POINTER is type, lowercase is object
lib.czeta.restype = c_void_p
lib.set_dgam.restype = None
lib.set_lm.restype = None
lib.zeta_eval.restype = ctypes.POINTER(c_double)

class czeta(object):
  def __init__(self):
    self.obj = c_void_p( lib.czeta())

  def set_dgam(self,dx,dy,dz,gam):
    print gam
    print c_double(gam)
    lib.set_dgam( self.obj, c_int(dx), c_int(dy), c_int(dz), c_double(gam))
    print "done"

  def set_lm(self,l,m):
    lib.set_lm(self.obj, c_int(l), c_int(m))

  def evaluate(self,q2):
    dout = lib.zeta_eval(self.obj,c_double(q2))
    return complex(dout[0],dout[1])

