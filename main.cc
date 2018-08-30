#include <cstdlib>
#include <exception>
#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_result.h>

// need gsl, blas => sudo apt-get install libgsl-dev libblas-dev

using namespace std;

struct gsl_underflow_exception : public exception {
   const char * what () const throw () {
      return "GSL_UNDERFLOW exception";
   }
};

struct gsl_other_exception : public exception {
   const char * what () const throw () {
      return "GSL_* exception";
   }
};

struct z00_params { double nx; double ny; double nz; double q2; };
double z00 (double x, void * p)
{
  gsl_sf_result res;
  struct z00_params * params = (struct z00_params *)p;
  double nx = (params->nx);
  double ny = (params->ny);
  double nz = (params->nz);
  double q2 = (params->q2);
  double arg = q2 - (nx*nx +ny*ny +nz*nz);
  // underflows are okay, just give 0. good enough.
  try {
    return gsl_sf_exp(x*arg); // exp{t *(q.q - n.n)}
  }
  catch ( gsl_underflow_exception& e) {
    return 0.;
  }
}

void gsl_to_c_handler(const char* reason, const char* file, int line, int gsl_errno)
{
  // throw my own error to catch
  if (gsl_errno == 15) { // don't know where to find other error codes
    throw gsl_underflow_exception();
  }
  std::cout << "reason   : " << reason    << std::endl;
  std::cout << "file     : " << file      << std::endl;
  std::cout << "line     : " << line      << std::endl;
  std::cout << "gsl_errno: " << gsl_errno << std::endl;
  std::cout << "GSL error: " << gsl_strerror( gsl_errno) << std::endl;
  throw gsl_other_exception();
}

int main(int argc, char** argv)
{

  double epsabs = 1e-8;
  double epsrel = 0.;
  double abserr = 0.;
  double result = 0.;
  size_t limit = 100;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(100);
  struct z00_params params = { 1.0, 0.0, 0.0, 0.0 };
  F.function = &z00;
  F.params = &params;

  gsl_set_error_handler( &gsl_to_c_handler );
  //std::cout << "test: " << z00( -1e10, &params) << std::endl; // overflow, test
  //std::cout << "test: " << z00( 1e10, &params) << std::endl; // underflow, test
  //std::cout << "test: " << z00( 10., &params) << std::endl;
  gsl_integration_qagiu(&F, 1., epsabs, epsrel, limit, w, &result, &abserr);
  gsl_set_error_handler( NULL );

  std::cout << "result: " << result << std::endl;
  std::cout << "test successful" << std::endl;

  //if (m % 2 == 0)
  //expected = M_SQRTPI + gsl_sf_gamma(0.5*(1.0 + m));
  //else
  //expected = M_SQRTPI;
  //printf ("m = %d\n", m);
  //printf ("intervals = %zu\n", gsl_integration_fixed_n(w));
  //printf ("result = % .18f\n", result);
  //printf ("exact result = % .18f\n", expected);
  //printf ("actual error = % .18f\n", result - expected);
  //gsl_integration_fixed_free (w);

  return 0;
}
