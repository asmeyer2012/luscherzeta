#include <assert.h>
//#include <math.h>
#include <gsl/gsl_complex.h>      // definitions for complex numbers
#include <gsl/gsl_complex_math.h> // math with complex numbers
//#include <gsl/gsl_integration.h> // GSL integration
//#include <gsl/gsl_math.h>        // basic GSL math
//#include <gsl/gsl_sf_dawson.h>   // D(x) = \sqrt{\pi}{2} e^{-x^2} Erfi(x)
//#include <gsl/gsl_sf_exp.h>      // GSL exponential function

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "gsl_exception_handler.h"
#include "spherical.h"
#include "zeta.h"

// C++ definition of zeta function
class cxxzeta {
 public:
 cxxzeta();
 void set_svec_gamma(double sx, double sy, double sz, double gamma);
 void set_lm(int l, int m);
 void evaluate(double q2, double& reZeta, double& imZeta);
 private:
 spherical_harmonic sharm;
 struct full_params p;
};

cxxzeta::cxxzeta() {
 this->sharm = spherical_harmonic();
 this->p.sharm = &(this->sharm);
 this->p.sx = 0.;
 this->p.sy = 0.;
 this->p.sz = 0.;
 this->p.gamma = 1.0;
 this->p.l = 0;
 this->p.m = 0;
 this->p.u2 = 1e-1;
}

void cxxzeta::set_svec_gamma(double sx, double sy, double sz, double gamma) {
 this->p.sx = sx;
 this->p.sy = sy;
 this->p.sz = sz;
 this->p.gamma = gamma;
}

void cxxzeta::set_lm(int l, int m) {
 this->p.l = l;
 this->p.m = m;
}

void cxxzeta::evaluate(double u2, double& reZeta, double& imZeta)
{
  this->p.u2 = u2;
  gsl_set_error_handler( &gsl_to_c_handler ); // set my own error handler
  gsl_complex result = full_zeta_lm( this->p );
  reZeta = GSL_REAL(result);
  imZeta = GSL_IMAG(result);
  gsl_set_error_handler( NULL );
}

// provide external C definitions, callable by python with ctypes
extern "C" {
    void* czeta() { return ((void*) new cxxzeta()); }
    void set_svec_gamma(void* cz, double sx, double sy, double sz, double gamma) {
      ((cxxzeta*)cz)->set_svec_gamma( sx, sy, sz, gamma); }
    void set_lm(void* cz, int l, int m) {
      ((cxxzeta*)cz)->set_lm( l,m); }
    //void zeta_eval(cxxzeta* cz, double q2, double& reZeta, double& imZeta) {
    //  cz->evaluate( q2, reZeta, imZeta); }
    bool zeta_eval(void* cz, double u2, double *out) {
      double reZeta, imZeta;
      try {
        ((cxxzeta*)cz)->evaluate( u2, reZeta, imZeta);
        out[0] = reZeta; out[1] = imZeta;
        return true;
      } catch ( gsl_underflow_exception& e ) {
        out[0] = 0.; out[1] = 0.;
        return true;
      } catch ( gsl_other_exception& e ) {
        out[0] = 0.; out[1] = 0.;
        return false;
      }
    }
}

