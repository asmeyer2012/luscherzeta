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

  return 0;
}
