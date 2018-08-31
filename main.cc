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

// full parameters to get all prefactors correct
struct full_params { int dx; int dy; int dz; int nx; int ny; int nz; double q2; double gam; };
// just compute the parameters that are needed
struct zeta_params { double q2; double ngam2; };

struct zeta_params full_to_zeta_params( const struct full_params in_params ){
  struct zeta_params out_params;
  out_params.q2 = in_params.q2;
  double nx = double(in_params.nx);
  double ny = double(in_params.ny);
  double nz = double(in_params.nz);
  double n2 = (nx*nx +ny*ny +nz*nz);
  if (in_params.dx || in_params.dy || in_params.dz) {
    double dx = double(in_params.dx);
    double dy = double(in_params.dy);
    double dz = double(in_params.dz);
    // vector dot products
    double nd = (nx*dx +ny*dy +nz*dz);
    double d2 = (dx*dx +dy*dy +dz*dz);
    double gam = in_params.gam;
    // scale \vec{n} parallel to \vec{d} by \gamma, then square
    // include a factor of \pi^2 for convenience
    out_params.ngam2 = M_PI*M_PI*((gam*gam -1.) *(nd*nd/d2) +n2);
  } else {
    out_params.ngam2 = M_PI*M_PI*n2;
  }
  return out_params;
}

double integral_zeta_00 (double x, void * p)
{
  gsl_sf_result res;
  struct zeta_params * params = (struct zeta_params *)p;
  double q2 = (params->q2);
  double ngam2 = (params->ngam2);
  // underflows are okay, just give 0. good enough.
  try {
    return pow( M_PI/x, 1.5) *gsl_sf_exp(q2*x -ngam2/x); // exp{t *(q.q - n.n)}
  }
  catch ( gsl_underflow_exception& e) {
    return 0.;
  }
}

int main(int argc, char** argv)
{

  gsl_set_error_handler( &gsl_to_c_handler ); // set my own error handler

  struct full_params fparams;
  struct zeta_params zparams;
  fparams.dx = 0;
  fparams.dy = 0;
  fparams.dz = 0;
  fparams.nx = 1;
  fparams.ny = 0;
  fparams.nz = 0;
  fparams.q2 = 0.;
  fparams.gam = 1.;
  zparams = full_to_zeta_params( fparams);

  gsl_function F;
  F.function = &integral_zeta_00;
  F.params = &zparams;

  double epsabs = 0.;
  double epsrel = 1e-8;
  double abserr = 0.;
  double result = 0.;
  size_t limit = 100;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);
  gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);

  std::cout << "result: " << result << std::endl;
  std::cout << "error : " << abserr << std::endl;
  std::cout << "test successful" << std::endl;

  gsl_set_error_handler( NULL );
  return 0;
}
