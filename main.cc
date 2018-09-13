//#include <algorithm> // std::count
//#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "gsl_exception_handler.h"
#include "spherical.h"
#include "vector_sorting.h"
#include "zeta.h"

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

  gsl_set_error_handler( &gsl_to_c_handler ); // set my own error handler

  spherical_harmonic sharm(4);
  struct full_params fparams;
  struct zeta_params zparams;
  fparams.dx = 0;
  fparams.dy = 1;
  fparams.dz = 1;
  fparams.q2 = .9;
  fparams.gam = 1.1;
  fparams.sharm = &sharm;
  fparams.l = 0;
  fparams.m = 0;

  std::cout << "zeta 00 : " << full_zeta_med_00 (fparams) << std::endl;
  std::cout << "zeta 00 : " << full_zeta_00 (fparams) << std::endl;
  std::cout << "zeta 00 : " << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
  //for (int l = 0; l < 4; l++) {
  //  fparams.l = l;
  //  for (int m = -l; m < l+1; m++) {
  //    fparams.m = m;
  //    gsl_complex result = full_zeta_lm(fparams);
  //    std::cout << "zeta " <<l <<", " <<m <<": " << gsl_complex_to_string(result) << std::endl;
  //      //<< gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
  //  }
  //}

  fparams.l = 1;
  fparams.m = 0;
  std::cout << "zeta 1, 0 : " << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
  fparams.m = 1;
  std::cout << "zeta 1,+1 : " << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
  fparams.m = -1;
  std::cout << "zeta 1,-1 : " << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;

  //fparams.l = 2;
  //fparams.m = 0;
  //std::cout << "zeta 2, 0 : " << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
  //fparams.m = 1;
  //std::cout << "zeta 2,+1 : " << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
  //fparams.m = -1;
  //std::cout << "zeta 2,-1 : " << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
  //fparams.m = 2;
  //std::cout << "zeta 2,+2 : " << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;
  //fparams.m = -2;
  //std::cout << "zeta 2,-2 : " << gsl_complex_to_string(full_zeta_lm (fparams)) << std::endl;

  //associated_legendre algen(); // generator for associated legendre polynomials

  gsl_set_error_handler( NULL );
  return 0;
}
