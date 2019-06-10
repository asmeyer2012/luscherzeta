import ctypes
from ctypes import c_bool,c_int,c_double,c_void_p
#import numpy as np

lib = ctypes.cdll.LoadLibrary('./zeta/libczeta.so')
#lib = ctypes.cdll.LoadLibrary('/sdcc/u/ameyer/code/luscherzeta/zeta/libczeta.so')

lib.czeta.argtypes = [ ]
lib.set_dgam.argtypes = [ c_void_p, c_int, c_int, c_int, c_double]
lib.set_lm.argtypes = [ c_void_p, c_int, c_int]
lib.zeta_eval.argtypes = [ c_void_p, c_double, ctypes.POINTER(c_double)]
## uppercase POINTER is type, lowercase is object
lib.czeta.restype = c_void_p
lib.set_dgam.restype = None
lib.set_lm.restype = None
lib.zeta_eval.restype = c_bool

class czeta(object):
  def __init__(self):
    self.obj = c_void_p( lib.czeta())
    self.mem = {} ## save previous values and reuse solutions for speed
    self.cur_l = (0,)
    self.cur_m = (0,)
    self.cur_d = (0,0,0)
    self.cur_g = ('1.0',)
    self.cur_q2 = ('1.0',)

  def set_dgam(self,dx,dy,dz,gam):
    lib.set_dgam( self.obj, c_int(dx), c_int(dy), c_int(dz), c_double(gam))
    self.cur_d = (dx,dy,dz)
    self.cur_g = (str(gam),)

  def set_lm(self,l,m):
    lib.set_lm(self.obj, c_int(l), c_int(m))
    self.cur_l = (l,)
    self.cur_m = (m,)

  def evaluate(self,q2):
    key = self.cur_d +self.cur_g +self.cur_l +self.cur_m +(str(q2),)
    if key in self.mem.keys():
      return self.mem[ key]
    dout = (c_double*2)()
    bout = lib.zeta_eval(self.obj,c_double(q2),dout)
    if bout:
      self.mem[ key] = complex(dout[0],dout[1])
      return complex(dout[0],dout[1])
    else:
      raise ValueError("GSL Error")

