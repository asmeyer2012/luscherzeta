#ifndef __ZETA_WRAPPER_H__
#define __ZETA_WRAPPER_H__

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

class czeta {
 public:
 czeta();
 void set_dgam(int dx, int dy, int dz, double gam);
 void set_lm(int l, int m);
 void evaluate(double q2, double& reZeta, double& imZeta);
 private:
 spherical_harmonic sharm;
 struct full_params p;
};

czeta::czeta() {
 this->sharm = spherical_harmonic();
 this->p.sharm = &(this->sharm);
 this->p.dx = 0;
 this->p.dy = 1;
 this->p.dz = 1;
 this->p.gam = 1.1;
 this->p.l = 0;
 this->p.m = 0;
 this->p.q2 = 1e-1;
}

void czeta::set_dgam(int dx, int dy, int dz, double gam) {
 this->p.dx = dx;
 this->p.dy = dy;
 this->p.dz = dz;
 this->p.gam = gam;
}

void czeta::set_lm(int l, int m) {
 this->p.l = l;
 this->p.m = m;
}

void czeta::evaluate(double q2, double& reZeta, double& imZeta)
{
  this->p.q2 = q2;
  gsl_set_error_handler( &gsl_to_c_handler ); // set my own error handler
  gsl_complex result = full_zeta_lm( this->p );
  reZeta = GSL_REAL(result);
  imZeta = GSL_IMAG(result);
  gsl_set_error_handler( NULL );
}

#endif
