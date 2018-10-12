//#include <algorithm> // std::count
//#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "gsl_exception_handler.h"
#include "spherical.h"
#include "zeta.h"
#include "zeta_wrapper.h"

// need gsl, blas => sudo apt-get install libgsl-dev libblas-dev

using namespace std;

std::string gsl_complex_to_string( gsl_complex val )
{
  stringstream sout;
  sout <<"(" <<GSL_REAL(val) <<", " <<GSL_IMAG(val) <<")";
  return sout.str();
}

int main(int argc, char** argv)
{
  int maxl = 4;

  spherical_harmonic sharm(maxl);
  struct full_params fparams;
  fparams.dx = 0;
  fparams.dy = 1;
  fparams.dz = 1;
  fparams.q2 = .9;
  fparams.gam = 1.1;
  fparams.sharm = &sharm;
  fparams.l = 0;
  fparams.m = 0;

  czeta cz;
  cz.set_dgam( fparams.dx, fparams.dy, fparams.dz, fparams.gam);
  double reZeta,imZeta;

  for (int l = 0; l < maxl+1; l++) {
    for (int m = -l; m < l+1; m++) {
      fparams.l = l;
      fparams.m = m;
      cz.set_lm( fparams.l, fparams.m);
      cz.evaluate( fparams.q2, reZeta, imZeta );

      std::cout << "zeta " <<l <<", " <<m <<": "
        << reZeta <<", " <<imZeta << std::endl;

      //gsl_set_error_handler( &gsl_to_c_handler ); // set my own error handler
      //std::cout << "zeta " <<l <<", " <<m <<": "
      //  << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
      //gsl_set_error_handler( NULL );
    }
  }

  fparams.dx = 0;
  fparams.dy = 0;
  fparams.dz = 0;
  fparams.q2 = .9;
  fparams.gam = 1.0;
  cz.set_dgam( fparams.dx, fparams.dy, fparams.dz, fparams.gam);

  // test some specific ratios of zeta functions
  fparams.l = 4;
  fparams.m = 0;
  cz.set_lm( fparams.l, fparams.m);
  cz.evaluate( fparams.q2, reZeta, imZeta );
  double z40 = reZeta;
  fparams.l = 4;
  fparams.m = 4;
  cz.set_lm( fparams.l, fparams.m);
  cz.evaluate( fparams.q2, reZeta, imZeta );
  double z44 = reZeta;
  std::cout << "test ratio (z44/z40): " <<z44/z40 <<"," <<sqrt(70.)/14. <<std::endl;

  fparams.l = 6;
  fparams.m = 0;
  cz.set_lm( fparams.l, fparams.m);
  cz.evaluate( fparams.q2, reZeta, imZeta );
  double z60 = reZeta;
  fparams.l = 6;
  fparams.m = 4;
  cz.set_lm( fparams.l, fparams.m);
  cz.evaluate( fparams.q2, reZeta, imZeta );
  double z64 = reZeta;
  std::cout << "test ratio (z64/z60): " <<z64/z60 <<"," <<-sqrt(14.)/2. <<std::endl;

  fparams.l = 8;
  fparams.m = 0;
  cz.set_lm( fparams.l, fparams.m);
  cz.evaluate( fparams.q2, reZeta, imZeta );
  double z80 = reZeta;
  fparams.l = 8;
  fparams.m = 4;
  cz.set_lm( fparams.l, fparams.m);
  cz.evaluate( fparams.q2, reZeta, imZeta );
  double z84 = reZeta;
  std::cout << "test ratio (z84/z80): " <<z84/z80 <<"," <<sqrt(154.)/33. <<std::endl;

  fparams.l = 8;
  fparams.m = 8;
  cz.set_lm( fparams.l, fparams.m);
  cz.evaluate( fparams.q2, reZeta, imZeta );
  double z88 = reZeta;
  std::cout << "test ratio (z88/z80): " <<z88/z80 <<"," <<sqrt(1430.)/66. <<std::endl;

  return 0;
}
