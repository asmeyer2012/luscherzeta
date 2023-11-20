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

int fill_weight_int( int zero_weight, int Nrmax){
  if (rand() % 100 < zero_weight) {
    return 0.;
  } else {
    return (rand() % Nrmax);
  }
}

double fill_weight( int zero_weight, int Nrmax){
  if (rand() % 100 < zero_weight) {
    return 0.;
  } else {
    return double( (rand() % Nrmax) +1) /double( Nrmax);
  }
}

bool test_spherical_harmonics()
{
  int maxl = 4;

  spherical_harmonic sharm(maxl);

  int Nrand = 10000;
  int Nrmax = 10000;
  int zero_weight = 10;
  std::vector<double> x0, x1, x2;
  for (int i=0; i<Nrand; i++) {
    x0.push_back( fill_weight( zero_weight, Nrmax));
    x1.push_back( fill_weight( zero_weight, Nrmax));
    x2.push_back( fill_weight( zero_weight, Nrmax));
  }
  for (int i=0; i<Nrand; i++) {
    double x0_i = x0[ i];
    double x1_i = x1[ i];
    double x2_i = x2[ i];
    // test identity Y_{lm}(theta, phi) = (-1)^{m} Y*_{lm}(theta, phi)
    for (int l=0; l<maxl; l++) {
      for (int m=-l; m<l+1; m++) {
        double phase = ( (m % 2 == 0) ? 1. : -1. );
        gsl_complex eval0 = sharm.evaluate( l, m, x0_i, x1_i, x2_i);
        gsl_complex eval1 = gsl_complex_mul_real(
          gsl_complex_conjugate( sharm.evaluate( l, -m, x0_i, x1_i, x2_i)), phase);
        if (!is_close( eval0, eval1)) {
          std::cout << "spherical harmonics failed: (l,m) = (" <<l <<", " <<m <<"): "
          <<gsl_complex_to_string( eval0) <<", " <<gsl_complex_to_string( eval1) <<std::endl;
          return false;
        }
        if (m == 0 && !is_small( GSL_IMAG( eval0))) {
          std::cout << "spherical harmonics m=0 failed: (l,m) = (" <<l <<", " <<m <<"): "
          <<gsl_complex_to_string( eval0) <<std::endl;
          return false;
        }
      }
    }
  }
  return true;
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

bool test_zeta_params()
{
  for (int i=0; i<21; i++) { rand(); } // preincrement the RNG
  int maxl = 4;
  int Nrand = 100;
  int Nvrand = 30;
  //int Nrmax = 10000;
  int zero_weight = 0;
  double L_a = 48.;
  // it seems this is (mildly) sensitive to whether the parameters are sensible or not
  // to avoid these kinds of issues, try to generate some reasonable parameters
  std::vector< std::vector< int>> n0, n1, n2;
  //std::vector< std::vector< int>> n1;
  //std::vector< std::vector< int>> n2;
  std::vector<double> s0, s1, s2;
  std::vector<double> gamma, u2vals;
  for (int i=0; i<Nrand; i++) {
    n0.push_back( {});
    n1.push_back( {});
    n2.push_back( {});
    for (int j=0; j<Nvrand; j++) {
      int n0_i = fill_weight_int( zero_weight, 6);
      int n1_i = fill_weight_int( zero_weight, 6);
      int n2_i = fill_weight_int( zero_weight, 6);
      n0[i].push_back( n0_i);
      n1[i].push_back( n1_i);
      n2[i].push_back( n2_i);
    }
    int d0_i = fill_weight_int( zero_weight, 3);
    int d1_i = fill_weight_int( zero_weight, 3);
    int d2_i = fill_weight_int( zero_weight, 3);
    //int d0_i = 0;
    //int d1_i = 0;
    //int d2_i = 1;
    double aM1_i = ( double( rand()) /double( RAND_MAX));
    double aM2_i = aM1_i *(0.3 + 0.7* double( rand()) /double( RAND_MAX));
    double dsq = double( d0_i*d0_i +d1_i*d1_i +d2_i*d2_i);
    double aP2lab_i = (2. *M_PI /L_a) *(2. *M_PI /L_a) *dsq;
    double aEcom_i = (aM1_i +aM2_i) *(1.01 +0.1 *( double( rand()) /double( RAND_MAX)));
    double aElab_i = sqrt( aEcom_i *aEcom_i +aP2lab_i);
    //double aElab_i = (aM1_i +aM2_i) *(1.01
    //  +(2. *M_PI /L_a) *sqrt( dsq) *(1.0 +0.5 *double( rand()) /double( RAND_MAX)));
    //double aEcom_i = sqrt( aElab_i *aElab_i -aP2lab_i);
    double aE2com_i = aEcom_i *aEcom_i;
    double d_fac = (1. +(aM1_i *aM1_i -aM2_i *aM2_i) /aE2com_i);
    double s0_i = d_fac *d0_i;
    double s1_i = d_fac *d1_i;
    double s2_i = d_fac *d2_i;
    double gamma_i = aElab_i /sqrt( aE2com_i);
    double u2_i = (0.5 *L_a /M_PI) *(0.5 *L_a /M_PI) *(
     .25 *aE2com_i -.5 *(aM1_i *aM1_i +aM2_i *aM2_i)
     +.25 *(aM1_i *aM1_i -aM2_i *aM2_i) *(aM1_i *aM1_i -aM2_i *aM2_i) /aE2com_i );
    s0.push_back( s0_i);
    s1.push_back( s1_i);
    s2.push_back( s2_i);
    gamma.push_back( gamma_i);
    u2vals.push_back( u2_i);
    //std::cout <<"d " <<d0_i <<"," <<d1_i <<"," <<d2_i <<std::endl;
    //std::cout <<"aM1 " <<aM1_i <<", aM2 " <<aM2_i <<", aElab " <<aElab_i <<std::endl;
    //std::cout <<"aP2lab_i " <<aP2lab_i <<", aE2com_i " <<aE2com_i <<std::endl;
    //std::cout <<"s " <<s0_i <<"," <<s1_i <<"," <<s2_i <<std::endl;
    //std::cout <<"gamma " <<gamma_i <<", u2 " <<u2_i <<std::endl;
  }

  bool test_result = true;
  double res_real0, res_imag0;
  double res_real1, res_imag1;
  // for comparisons
  double epsilon_rel = 1e-8;
  double epsilon_abs = 4e-9;
  // for integration
  double epsabs = 0.;
  double epsrel = 1e-7;
  size_t limit = 1000;
  // integration error
  double abserr0 = 0;
  double abserr1 = 0;
  spherical_harmonic sharm( maxl);
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);
  F.function = &integral_zeta_lm;
  for (int i=0; i<Nrand; i++) {
    full_params fp0, fp1;
    fp0.sx = fp1.sx = s0[i];
    fp0.sy = fp1.sy = s0[i];
    fp0.sz = fp1.sz = s0[i];
    fp0.u2 = fp1.u2 = u2vals[i];
    fp0.gamma = fp1.gamma = gamma[i];
    fp0.sharm = fp1.sharm = &sharm;
    for (int j=0; j<Nvrand; j++) {
      for (int l=0; l<maxl; l++) {
        for (int m=0; m<l+1; m++) {
          fp0.l = fp1.l = l;
          fp0.m =  m;
          fp1.m = -m;
          zeta_lm_params fz0 = full_to_zeta_lm_params( n0[i][j], n1[i][j], n2[i][j], fp0);
          zeta_lm_params fz1 = full_to_zeta_lm_params( n0[i][j], n1[i][j], n2[i][j], fp1);

          res_real0 = GSL_REAL( fz0.rlYlm);
          res_imag0 = GSL_IMAG( fz0.rlYlm);
          res_real1 = GSL_REAL( fz1.rlYlm);
          res_imag1 = GSL_IMAG( fz1.rlYlm);
          // Yl-m = (-1)^m Ylm*
          if ((m % 2) == 0) { res_imag1 *= -1; }
          else              { res_real1 *= -1; }
          if (!is_close( res_real0, res_real1, epsilon_rel, epsilon_abs)
              || !is_close( res_imag0, res_imag1, epsilon_rel, epsilon_abs)) {
            std::cout <<"s " <<s0[i] <<"," <<s1[i] <<"," <<s2[i] <<std::endl;
            std::cout <<"gamma " <<gamma[i] <<", u2 " <<u2vals[i] <<std::endl;
            std::cout <<"n " <<n0[i][j] <<"," <<n1[i][j] <<"," <<n2[i][j] <<std::endl;
            std::cout << "zeta_params parity rlYlm failed: (l,m) = (" <<l <<", " <<m <<"): ("
              <<res_real0 <<", " <<res_imag0 <<"), ("
              <<res_real1 <<", " <<res_imag1 <<")" <<std::endl;
            test_result = false;
          }

          res_real0 = GSL_REAL( fz0.wlYlm);
          res_imag0 = GSL_IMAG( fz0.wlYlm);
          res_real1 = GSL_REAL( fz1.wlYlm);
          res_imag1 = GSL_IMAG( fz1.wlYlm);
          // Yl-m = (-1)^m Ylm*
          if ((m % 2) == 0) { res_imag1 *= -1; }
          else              { res_real1 *= -1; }
          // counteract the -i factor in wlYlm
          if ((l % 2) == 1) { res_real1 *= -1; res_imag1 *= -1; }
          if (!is_close( res_real0, res_real1, epsilon_rel, epsilon_abs)
              || !is_close( res_imag0, res_imag1, epsilon_rel, epsilon_abs)) {
            std::cout <<"s " <<s0[i] <<"," <<s1[i] <<"," <<s2[i] <<std::endl;
            std::cout <<"gamma " <<gamma[i] <<", u2 " <<u2vals[i] <<std::endl;
            std::cout <<"n " <<n0[i][j] <<"," <<n1[i][j] <<"," <<n2[i][j] <<std::endl;
            std::cout << "zeta_params parity wlYlm failed: (l,m) = (" <<l <<", " <<m <<"): ("
              <<res_real0 <<", " <<res_imag0 <<"), ("
              <<res_real1 <<", " <<res_imag1 <<")" <<std::endl;
            test_result = false;
          }

          res_real0 = GSL_REAL( fz0.leadsum);
          res_imag0 = GSL_IMAG( fz0.leadsum);
          res_real1 = GSL_REAL( fz1.leadsum);
          res_imag1 = GSL_IMAG( fz1.leadsum);
          // Yl-m = (-1)^m Ylm*
          if ((m % 2) == 0) { res_imag1 *= -1; }
          else              { res_real1 *= -1; }
          if (!is_close( res_real0, res_real1, epsilon_rel, epsilon_abs)
              || !is_close( res_imag0, res_imag1, epsilon_rel, epsilon_abs)) {
            std::cout <<"s " <<s0[i] <<"," <<s1[i] <<"," <<s2[i] <<std::endl;
            std::cout <<"gamma " <<gamma[i] <<", u2 " <<u2vals[i] <<std::endl;
            std::cout <<"n " <<n0[i][j] <<"," <<n1[i][j] <<"," <<n2[i][j] <<std::endl;
            std::cout << "zeta_params parity leadsum failed: (l,m) = (" <<l <<", " <<m <<"): ("
              <<res_real0 <<", " <<res_imag0 <<"), ("
              <<res_real1 <<", " <<res_imag1 <<")" <<std::endl;
            test_result = false;
          }

          gsl_set_error_handler( &gsl_to_c_handler );
          // skip n=(0,0,0) since it produces infinity
          bool integration_success = !(n0[i][j] == 0 && n1[i][j] == 0 && n2[i][j] == 0);
          try {
            if (integration_success) {
              F.params = &fz0;
              gsl_integration_qag( &F, 0., 1., epsabs, epsrel, limit, 1, w, &res_real0, &abserr0);
              integration_success = true;
            }
          }
          catch ( gsl_other_exception& e ) {
            std::cout <<"s " <<s0[i] <<"," <<s1[i] <<"," <<s2[i] <<std::endl;
            std::cout <<"gamma " <<gamma[i] <<", u2 " <<u2vals[i] <<std::endl;
            std::cout <<"n " <<n0[i][j] <<"," <<n1[i][j] <<"," <<n2[i][j] <<std::endl;
            std::cout << "zeta_params parity integration failed: (l,m) = (" <<l <<", " <<m <<"): ("
              <<res_real0 <<")" <<std::endl;
            integration_success = false;
          }
          try {
            if (integration_success) {
              F.params = &fz1;
              gsl_integration_qag( &F, 0., 1., epsabs, epsrel, limit, 1, w, &res_real1, &abserr1);
            }
          }
          catch ( gsl_other_exception& e ) {
            std::cout <<"s " <<s0[i] <<"," <<s1[i] <<"," <<s2[i] <<std::endl;
            std::cout <<"gamma " <<gamma[i] <<", u2 " <<u2vals[i] <<std::endl;
            std::cout <<"n " <<n0[i][j] <<"," <<n1[i][j] <<"," <<n2[i][j] <<std::endl;
            std::cout << "zeta_params parity integration failed: (l,m) = (" <<l <<", " <<m <<"): ("
              <<res_real1 <<")" <<std::endl;
          }
          gsl_set_error_handler( NULL );
          if (integration_success) {
            if (!is_close( res_real0, res_real1, epsilon_rel, epsilon_abs)) {
              std::cout <<"s " <<s0[i] <<"," <<s1[i] <<"," <<s2[i] <<std::endl;
              std::cout <<"gamma " <<gamma[i] <<", u2 " <<u2vals[i] <<std::endl;
              std::cout <<"n " <<n0[i][j] <<"," <<n1[i][j] <<"," <<n2[i][j] <<std::endl;
              std::cout << "zeta_params integration check failed: (l,m) = (" <<l <<", " <<m <<"): ("
                <<res_real0 <<"), (" <<res_real1 <<")" <<std::endl;
              test_result = false;
            }
          }

          // for iphase, sensitive to Ylm * sin(n.s)
          // take n => -n to do proper comparison
          fz1 = full_to_zeta_lm_params( -n0[i][j], -n1[i][j], -n2[i][j], fp1);
          res_real0 = GSL_REAL( fz0.iphase);
          res_imag0 = GSL_IMAG( fz0.iphase);
          res_real1 = GSL_REAL( fz1.iphase);
          res_imag1 = GSL_IMAG( fz1.iphase);
          // Yl-m = (-1)^m Ylm*
          if ((m % 2) == 0) { res_imag1 *= -1; }
          else              { res_real1 *= -1; }
          //// spherical harmonics also satisfy Ylm(-r) = (-1)^l Ylm(r)
          //if ((l % 2) == 1) { res_real1 *= -1; res_imag1 *= -1; }
          if (!is_close( res_real0, res_real1, epsilon_rel, epsilon_abs)
              || !is_close( res_imag0, res_imag1, epsilon_rel, epsilon_abs)) {
            std::cout <<"s " <<s0[i] <<"," <<s1[i] <<"," <<s2[i] <<std::endl;
            std::cout <<"gamma " <<gamma[i] <<", u2 " <<u2vals[i] <<std::endl;
            std::cout <<"n " <<n0[i][j] <<"," <<n1[i][j] <<"," <<n2[i][j] <<std::endl;
            std::cout << "zeta_params parity iphase failed: (l,m) = (" <<l <<", " <<m <<"): ("
              <<res_real0 <<", " <<res_imag0 <<"), ("
              <<res_real1 <<", " <<res_imag1 <<")" <<std::endl;
            test_result = false;
          }

        }
      }
    }
  }
  return test_result;
}

bool test_zeta_parity( czeta &cz)
{
  for (int i=0; i<29; i++) { rand(); } // preincrement the RNG
  int maxl = 4;
  int Nrand = 20;
  //int Nrmax = 10000;
  int zero_weight = 0;
  double L_a = 48.;
  // try to generate some reasonable parameters
  std::vector<double> s0, s1, s2;
  std::vector<double> gamma, u2vals;
  for (int i=0; i<Nrand; i++) {
    int d0_i = fill_weight_int( zero_weight, 3);
    int d1_i = fill_weight_int( zero_weight, 3);
    int d2_i = fill_weight_int( zero_weight, 3);
    double aM1_i = ( double( rand()) /double( RAND_MAX));
    double aM2_i = aM1_i *(0.3 + 0.7* double( rand()) /double( RAND_MAX));
    double dsq = double( d0_i*d0_i +d1_i*d1_i +d2_i*d2_i);
    double aP2lab_i = (2. *M_PI /L_a) *(2. *M_PI /L_a) *dsq;
    double aEcom_i = (aM1_i +aM2_i) *(1.01 +0.1 *( double( rand()) /double( RAND_MAX)));
    double aElab_i = sqrt( aEcom_i *aEcom_i +aP2lab_i);
    double aE2com_i = aEcom_i *aEcom_i;
    double d_fac = (1. +(aM1_i *aM1_i -aM2_i *aM2_i) /aE2com_i);
    double s0_i = d_fac *d0_i;
    double s1_i = d_fac *d1_i;
    double s2_i = d_fac *d2_i;
    double gamma_i = aElab_i /sqrt( aE2com_i);
    double u2_i = (0.5 *L_a /M_PI) *(0.5 *L_a /M_PI) *(
     .25 *aE2com_i -.5 *(aM1_i *aM1_i +aM2_i *aM2_i)
     +.25 *(aM1_i *aM1_i -aM2_i *aM2_i) *(aM1_i *aM1_i -aM2_i *aM2_i) /aE2com_i );
    s0.push_back( s0_i);
    s1.push_back( s1_i);
    s2.push_back( s2_i);
    gamma.push_back( gamma_i);
    u2vals.push_back( u2_i);
    //std::cout <<"d " <<d0_i <<"," <<d1_i <<"," <<d2_i <<std::endl;
    //std::cout <<"aM1 " <<aM1_i <<", aM2 " <<aM2_i <<", aElab " <<aElab_i <<std::endl;
    //std::cout <<"aP2lab_i " <<aP2lab_i <<", aE2com_i " <<aE2com_i <<std::endl;
    //std::cout <<"s " <<s0_i <<"," <<s1_i <<"," <<s2_i <<std::endl;
    //std::cout <<"gamma " <<gamma_i <<", u2 " <<u2_i <<std::endl;
  }

  bool test_result = true;
  for (int i=0; i<Nrand; i++) {
    std::cout <<" -- test " <<i+1 <<"/" <<Nrand <<std::endl;
    cz.set_svec_gamma( s0[i], s1[i], s2[i], gamma[i]);
    double res_real0 = 0., res_real1 = 0.;
    double res_imag0 = 0., res_imag1 = 0.;
    //std::cout <<"s " <<s0[i] <<"," <<s1[i] <<"," <<s2[i] <<std::endl;
    //std::cout <<"gamma " <<gamma[i] <<", u2 " <<u2vals[i] <<std::endl;
    // test identity Z_{lm}(s, gamma, u2) = (-1)^{m} Z*_{lm}(s, gamma, u2)
    for (int l=0; l<maxl; l++) {
      for (int m=0; m<l+1; m++) {
        cz.set_lm( l, m);
        cz.evaluate( u2vals[i], res_real0, res_imag0);
        cz.set_lm( l, -m);
        cz.evaluate( u2vals[i], res_real1, res_imag1);
        // apply (-1)^m phase and complex conjugation
        if ((m % 2) == 0) { res_imag1 *= -1; }
        else              { res_real1 *= -1; }
        if (!is_close( res_real0, res_real1, 1e-6, 1e-8)
            || !is_close( res_imag0, res_imag1, 1e-6, 1e-8)) {
          std::cout << "zeta parity failed: (l,m) = (" <<l <<", " <<m <<"): ("
            <<res_real0 <<", " <<res_imag0 <<"), ("
            <<res_real1 <<", " <<res_imag1 <<")" <<std::endl;
          test_result = false;
        }
        else {
          //std::cout << "zeta parity success: (l,m) = (" <<l <<", " <<m <<"): ("
          //  <<res_real0 <<", " <<res_imag0 <<"), ("
          //  <<res_real1 <<", " <<res_imag1 <<")" <<std::endl;
        }
      }
    }
  }
  return test_result;
}

int main(int argc, char** argv)
{
  czeta cz;

  // initialize czeta calculator
  double sx = 0.;
  double sy = 0.;
  double sz = 0.;
  double gamma = 1.0;
  cz.set_svec_gamma( sx, sy, sz, gamma);

  bool test_success = true;

  if (test_success) {
    std::cout <<"spherical harmonic tests" <<std::endl;
    if (test_success) { test_success = test_spherical_harmonics(); }
    if (test_success) {
      std::cout <<"spherical harmonic tests passed" <<std::endl;
    } else {
      std::cout <<"spherical harmonic tests failed" <<std::endl;
    }
  }

  if (test_success) {
    std::cout <<"zeta_params tests" <<std::endl;
    if (test_success) { test_success = test_zeta_params(); }
    if (test_success) {
      std::cout <<"zeta_params tests passed" <<std::endl;
    } else {
      std::cout <<"zeta_params tests failed" <<std::endl;
    }
  }

  if (test_success) {
    std::cout <<"ratio_test_com tests" <<std::endl;
    if (test_success) { test_success = ratio_test_com( cz, .9); }
    if (test_success) { test_success = ratio_test_com( cz, .1234); }
    if (test_success) {
      std::cout <<"ratio_test_com passed" <<std::endl;
    } else {
      std::cout <<"ratio_test_com failed" <<std::endl;
    }
  }

  if (test_success) {
    std::cout <<"zeta parity tests" <<std::endl;
    if (test_success) { test_success = test_zeta_parity( cz); }
    if (test_success) {
      std::cout <<"zeta parity test passed" <<std::endl;
    } else {
      std::cout <<"zeta parity test failed" <<std::endl;
    }
  }

  //int maxl = 4;
  //double reZeta, imZeta;
  //sx = 0.;
  //sy = 1.;
  //sz = 1.;
  //gamma = 1.1;
  //cz.set_svec_gamma( sx, sy, sz, gamma);

  //double u2 = .9;
  //for (int l = 0; l < maxl+1; l++) {
  //  for (int m = -l; m < l+1; m++) {
  //    cz.set_lm( l, m);
  //    cz.evaluate( u2, reZeta, imZeta );

  //    std::cout << "zeta " <<l <<", " <<m <<": "
  //      << reZeta <<", " <<imZeta << std::endl;

  //    //gsl_set_error_handler( &gsl_to_c_handler ); // set my own error handler
  //    //std::cout << "zeta " <<l <<", " <<m <<": "
  //    //  << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
  //    //gsl_set_error_handler( NULL );
  //  }
  //}

  return 0;
}
