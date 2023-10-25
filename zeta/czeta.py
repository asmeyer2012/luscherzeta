import ctypes
from ctypes import c_bool,c_int,c_double,c_void_p

lib = ctypes.cdll.LoadLibrary('./zeta/libczeta.so')

lib.czeta.argtypes = [ ]
lib.set_svec_gamma.argtypes = [ c_void_p, c_double, c_double, c_double, c_double]
lib.set_lm.argtypes = [ c_void_p, c_int, c_int]
lib.zeta_eval.argtypes = [ c_void_p, c_double, ctypes.POINTER( c_double)]
## uppercase POINTER is type, lowercase is object
lib.czeta.restype = c_void_p
lib.set_svec_gamma.restype = None
lib.set_lm.restype = None
lib.zeta_eval.restype = c_bool

class czeta(object):
  def __init__( self, do_caching=False):
    self.obj = c_void_p( lib.czeta())
    self.mem = {} ## save previous values and reuse solutions for speed
    self.cur_l = (0,)
    self.cur_m = (0,)
    self.cur_s = (0, 0, 0)
    self.cur_g = ('1.0',)
    self.cur_u2 = ('1.0',)
    self.do_caching = do_caching

  def set_svec_gamma( self, sx, sy, sz, gamma):
    lib.set_svec_gamma(
      self.obj,
      c_double( sx), c_double( sy), c_double( sz),
      c_double( gamma))
    self.cur_s = (sx, sy, sz)
    self.cur_g = (str(gamma),)

  def set_lm(self,l,m):
    lib.set_lm(self.obj, c_int(l), c_int(m))
    self.cur_l = (l,)
    self.cur_m = (m,)

  def evaluate( self, u2):
    if self.do_caching:
      key = self.cur_d +self.cur_g +self.cur_l +self.cur_m +(str(u2),)
      if key in self.mem.keys():
        return self.mem[ key]
    dout = (c_double*2)()
    bout = lib.zeta_eval( self.obj, c_double(u2), dout)
    if bout:
      if self.do_caching:
        self.mem[ key] = complex( dout[0], dout[1])
      return complex( dout[0], dout[1])
    else:
      raise ValueError("GSL Error")

