#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <float.h> // contains definitions of FLT_EPSILON and DBL_EPSILON

#include "gsl_exception_handler.h"
#include "zeta.h"
#include "zeta_wrapper.h"

// need gsl, blas => sudo apt-get install libgsl-dev libblas-dev

using namespace std;

/// Absolute size comparison against 0 for double
bool is_small( double x) { return fabs( x) < 16384 *DBL_EPSILON; }

/// Compare if two doubles are close in value to each other.
/**
  Will return true for numbers that share the first 10 significant digits
  or if both numbers are below about \f$10^{-12}\f$.
*/
bool is_close(
  double a, double b,
  double epsilon = 128 * DBL_EPSILON, double abs_th = 16384 *DBL_EPSILON)
{
  if (a == b) return true;

  auto diff = abs( a -b);
  auto norm = min( (abs(a) + abs(b)), std::numeric_limits<double>::max());
  // or even faster: std::min(std::abs(a + b), std::numeric_limits<float>::max());
  return diff < max( abs_th, epsilon * norm);
}

std::string gsl_complex_to_string( gsl_complex val )
{
  stringstream sout;
  sout <<"(" <<GSL_REAL(val) <<", " <<GSL_IMAG(val) <<")";
  return sout.str();
}

// code checks whether (real) eval of zeta function in center of mass gives expected ratios
// imaginary parts are assumed to be identically zero
bool ratio_test_com_single(
    czeta &cz, double u2, double ratio_val_re,
    int l0, int m0, int l1, int m1) {

  // initialize fparams
  cz.set_svec_gamma( 0., 0., 0., 1.);

  // for dumping results
  double reZeta0, imZeta0;
  double reZeta1, imZeta1;

  // test some specific ratios of zeta functions
  cz.set_lm( l0, m0);
  cz.evaluate( u2, reZeta0, imZeta0 );
  cz.set_lm( l1, m1);
  cz.evaluate( u2, reZeta1, imZeta1 );
  bool test_success = (
    is_close( reZeta0/reZeta1, ratio_val_re)
    && is_small( imZeta0) && is_small( imZeta1) );

  if (!test_success) {
    if (!is_small( imZeta0)) {
      std::cout << "imaginary test failed for eval 0: " <<imZeta0 <<"!=" <<0. <<std::endl;
    } else if (!is_small( imZeta1)) {
      std::cout << "imaginary test failed for eval 1: " <<imZeta1 <<"!=" <<0. <<std::endl;
    } else {
      std::cout << "test ratio failed: " <<reZeta0 /reZeta1 <<"!=" <<ratio_val_re <<std::endl;
    }
  }

  return test_success;
}

bool ratio_test_com( czeta &cz, double u2) {

  bool test_success = true;

  if (test_success) {
    test_success = ratio_test_com_single( cz, u2, sqrt(70.)/14., 4, 4, 4, 0);
    if (!test_success) {
      std::cout <<"z44/z40 test failed" <<std::endl;
    }
  }

  if (test_success) {
    test_success = ratio_test_com_single( cz, u2, -sqrt(14.)/2., 6, 4, 6, 0);
    if (!test_success) {
      std::cout <<"z64/z60 test failed" <<std::endl;
    }
  }

  if (test_success) {
    test_success = ratio_test_com_single( cz, u2, sqrt(154.)/33., 8, 4, 8, 0);
    if (!test_success) {
      std::cout <<"z84/z80 test failed" <<std::endl;
    }
  }

  if (test_success) {
    test_success = ratio_test_com_single( cz, u2, sqrt(1430.)/66., 8, 8, 8, 0);
    if (!test_success) {
      std::cout <<"z88/z80 test failed" <<std::endl;
    }
  }

  return test_success;
}

int main(int argc, char** argv)
{
  int maxl = 4;
  double reZeta, imZeta;
  czeta cz;

  // initialize czeta calculator
  double sx = 0.;
  double sy = 0.;
  double sz = 0.;
  double gamma = 1.0;
  cz.set_svec_gamma( sx, sy, sz, gamma);

  bool test_success = true;

  std::cout <<"ratio_test_com tests" <<std::endl;
  if (test_success) { test_success = ratio_test_com( cz, .9); }
  if (test_success) { test_success = ratio_test_com( cz, .1234); }
  if (test_success) {
    std::cout <<"ratio_test_com passed" <<std::endl;
  } else {
    std::cout <<"ratio_test_com failed" <<std::endl;
  }

  sx = 0.;
  sy = 1.;
  sz = 1.;
  gamma = 1.1;
  cz.set_svec_gamma( sx, sy, sz, gamma);

  double u2 = .9;
  for (int l = 0; l < maxl+1; l++) {
    for (int m = -l; m < l+1; m++) {
      cz.set_lm( l, m);
      cz.evaluate( u2, reZeta, imZeta );

      std::cout << "zeta " <<l <<", " <<m <<": "
        << reZeta <<", " <<imZeta << std::endl;

      //gsl_set_error_handler( &gsl_to_c_handler ); // set my own error handler
      //std::cout << "zeta " <<l <<", " <<m <<": "
      //  << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
      //gsl_set_error_handler( NULL );
    }
  }

  return 0;
}
