//#include <algorithm> // std::count
//#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <sstream>
//#include <vector>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_result.h>

#include "gsl_exception_handler.h"
#include "vector_sorting.h"

// need gsl, blas => sudo apt-get install libgsl-dev libblas-dev

using namespace std;

std::string gsl_complex_to_string( gsl_complex val )
{
  stringstream sout;
  sout <<"(" <<GSL_REAL(val) <<", " <<GSL_IMAG(val) <<")";
  return sout.str();
}

// full parameters to get all prefactors correct
struct full_params { int dx; int dy; int dz; int nx; int ny; int nz; double q2; double gam; };
// just compute the parameters that are needed
struct zeta_params {
 double q2;          // q2 from input
 double ngam2;       // exponent of integral
 double leadsum;     // leading sum term
 double iphase;      // prefactor and phase on integral
 };

// \sum_l (q^2)^l /((2l-1) \Gamma(l+1) ) == (\sqrt{q^2} Erfi{\sqrt{q^2}} - e^{q^2})
struct zeta_params full_to_zeta_params( const struct full_params in_params )
{
  struct zeta_params out_params;
  double nx = double(in_params.nx);
  double ny = double(in_params.ny);
  double nz = double(in_params.nz);
  double n2 = (nx*nx +ny*ny +nz*nz);
  out_params.q2 = in_params.q2;
  double r2q2;

  if (in_params.dx || in_params.dy || in_params.dz) {
    double gam = in_params.gam;
    double dx = double(in_params.dx);
    double dy = double(in_params.dy);
    double dz = double(in_params.dz);
    // vector dot products
    double nd = (nx*dx +ny*dy +nz*dz);
    double d2 = (dx*dx +dy*dy +dz*dz);
    // (n_par + d/2)^2/gam^2 +n_perp^2
    double r2 = (nd*nd/d2 +nd +.25*d2)/(gam*gam) +n2 - nd*nd/d2;
    std::cout << "nd :" << nd << std::endl;
    std::cout << "d2 :" << d2 << std::endl;
    std::cout << "r2 :" << r2 << std::endl;

    // scale \vec{n} parallel to \vec{d} by \gamma, then square
    // include a factor of \pi^2 for convenience
    out_params.ngam2 = M_PI*M_PI*((gam*gam -1.) *(nd*nd/d2) +n2);
    // == \gamma e^{-i\pi n.d} / 2\sqrt{\pi}, phase for integral
    out_params.iphase = .5*gam /M_SQRTPI *gsl_pow_int(-1.,int(nd)); // nd is always integer, so real
    r2q2 = r2 - in_params.q2; // argument for leading sum term
  }
  else {

    // \vec{d} == 0, \gamma == 1.
    out_params.ngam2 = M_PI*M_PI*n2;
    out_params.iphase = .5 /M_SQRTPI; // == 1. / 2\sqrt{\pi}
    r2q2 = n2 - in_params.q2; // == n^2 - q^2
  }

  out_params.leadsum = gsl_sf_exp(-r2q2) /(2.*M_SQRTPI*r2q2);
  std::cout << "ngam2  : " << out_params.ngam2  << std::endl;
  std::cout << "iphase : " << out_params.iphase  << std::endl;
  std::cout << "leadsum: " << out_params.leadsum << std::endl;
  return out_params;
}

double integral_zeta_00 (double x, void * p)
{
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
  fparams.dx = 1;
  fparams.dy = 1;
  fparams.dz = 0;
  fparams.nx = 3;
  fparams.ny = 0;
  fparams.nz = 0;
  fparams.q2 = .1;
  fparams.gam = 1.1;
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
  double cresult = zparams.iphase *result;

  std::cout << "result: " << cresult << std::endl;
  std::cout << "error : " << abserr  << std::endl;
  std::cout << "integral test successful" << std::endl;

  //std::vector< std::vector<int>> testvec;
  //testvec.push_back({0,0,1}); // done
  //testvec.push_back({0,1,1}); // done
  //testvec.push_back({1,1,1}); // done
  //testvec.push_back({0,1,2}); // done
  //testvec.push_back({1,1,2}); // done
  //testvec.push_back({1,2,2}); // done
  //testvec.push_back({1,2,3}); // done
  //for (auto vec = testvec.begin(); vec != testvec.end(); vec++) {
  //  auto vecout = all_permutations( *vec);
  //  std::cout << "input vector: " <<(*vec)[0] <<", " <<(*vec)[1] <<", " <<(*vec)[2] << std::endl;
  //  int i = 0;
  //  for (auto pvec = vecout.begin(); pvec != vecout.end(); pvec++) {
  //    i++;
  //    std::cout << "-- output vector " <<i <<": "
  //    <<(*pvec)[0] <<", " <<(*pvec)[1] <<", " <<(*pvec)[2] << std::endl;
  //  }
  //}

  for (int i=1; i<25; i++){
    auto vecout = all_combos( i);
    std::cout << "i: " <<i << std::endl;
    for (auto pvec = vecout.begin(); pvec != vecout.end(); pvec++) {
      std::cout << "-- output vector: "
      <<(*pvec)[0] <<", " <<(*pvec)[1] <<", " <<(*pvec)[2] << std::endl;
    }
  }

  gsl_set_error_handler( NULL );
  return 0;
}
