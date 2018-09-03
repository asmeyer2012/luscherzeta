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
#include <gsl/gsl_sf_dawson.h> // D(x) = \sqrt{\pi}{2} e^{-x^2} Erfi(x)
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_result.h>

#include "gsl_exception_handler.h"
#include "vector_sorting.h"

// need gsl, blas => sudo apt-get install libgsl-dev libblas-dev

#define RELERR 1e-8 // relative error to compute zeta functions to

using namespace std;

std::string gsl_complex_to_string( gsl_complex val )
{
  stringstream sout;
  sout <<"(" <<GSL_REAL(val) <<", " <<GSL_IMAG(val) <<")";
  return sout.str();
}

// complex error function
double erfi( double q2)
{ return 2. *gsl_sf_exp( gsl_pow_int( q2, 2)) *gsl_sf_dawson( q2) /M_SQRTPI; }

// check relative error of two values
bool relerr_check (double prevVal, double nextVal)
{ return ( RELERR *abs( nextVal) < abs( nextVal -prevVal) ); }

// full parameters to get all prefactors correct
struct full_params { int dx; int dy; int dz; double q2; double gam; };
// just compute the parameters that are needed
struct zeta_params {
 double r2q2;        // r2 - q2, for consistency checks
 double q2;          // q2 from input
 double ngam2;       // exponent of integral
 double leadsum;     // leading sum term
 double iphase;      // prefactor and phase on integral
 };

// \sum_l (q^2)^l /((2l-1) \Gamma(l+1) ) == (\sqrt{q^2} Erfi{\sqrt{q^2}} - e^{q^2})
struct zeta_params full_to_zeta_params( int inx, int iny, int inz,
  const struct full_params in_params )
{
  struct zeta_params out_params;
  double nx = double(inx);
  double ny = double(iny);
  double nz = double(inz);
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
    //std::cout << "nd :" << nd << std::endl;
    //std::cout << "d2 :" << d2 << std::endl;
    //std::cout << "r2 :" << r2 << std::endl;

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

  out_params.r2q2 = r2q2;
  out_params.leadsum = gsl_sf_exp(-r2q2) /(2.*M_SQRTPI*r2q2);
  //std::cout << "ngam2  : " << out_params.ngam2  << std::endl;
  //std::cout << "iphase : " << out_params.iphase  << std::endl;
  //std::cout << "leadsum: " << out_params.leadsum << std::endl;
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

// n = (0,0,0) divergent term, added separately from others
// total should be (1/\sqrt{4\pi}) 2\gamma \pi^{3/2} ( \sqrt{\pi q2} erfi(\sqrt{q2}) - exp{q2} )
double n000_sum( double q2)
{
  double srq2 = sqrt(q2);
  return M_SQRTPI *srq2 *erfi(srq2) -gsl_sf_exp(q2);
}

// very slowly converging!
double full_zeta_slow_00 (struct full_params p)
{
  int i = 1;
  double prevAns = 0.;
  double nextAns = 0.;
  bool skipP2 = false;
  // add the 000 term
  struct zeta_params zp = full_to_zeta_params( 0,0,0, p);
  nextAns += .5 /M_SQRTPI /zp.r2q2;
  //std::cout << "  |  contribution : " << 1./zp.r2q2 << std::endl;
  //while ( (RELERR*abs(nextAns) < abs(prevAns - nextAns)) || (i % 4 != 0) || skipP2 ) {
  while ( relerr_check( prevAns, nextAns) || (i % 4 != 0) || skipP2 ) {
    //std::cout << "iteration " <<i <<std::endl;
    prevAns = nextAns; skipP2 = false;
    auto vecCombos = all_combos( i); // get a list of all vector combos for this choice
    if ( vecCombos.size() == 0) { skipP2 = true; }
    for ( auto vecc = vecCombos.begin(); vecc != vecCombos.end(); vecc++ ) {
      auto vecPerms = all_permutations( *vecc);
      for ( auto vecp = vecPerms.begin(); vecp != vecPerms.end(); vecp++ ) {
        //std::cout << "sum vector: "
        //  <<(*vecp)[0] <<", " <<(*vecp)[1] <<", " <<(*vecp)[2];
        zp = full_to_zeta_params( (*vecp)[0], (*vecp)[1], (*vecp)[2], p);
        //nextAns += 1./zp.r2q2;
        nextAns += .5 /M_SQRTPI /zp.r2q2;
        std::cout << "  |  contribution : " << 1./zp.r2q2 << std::endl;
      }
    }
    //std::cout << "skipP2 : " << skipP2 << std::endl;
    //std::cout << "prevAns: " << prevAns << std::endl;
    //std::cout << "nextAns: " << nextAns << std::endl;
    //std::cout << RELERR << " < " << abs(prevAns - nextAns)/abs(nextAns) << std::endl;
    i++;
    assert(i < 100);
  }
  return nextAns;
}

// integration method, faster
// q2 = 0, gam = 1. => zeta_00 = -8.91363292
double full_zeta_00 (struct full_params p)
{
  int i = 1;
  double prevAnsSum = 0.;
  double nextAnsSum = 0.;
  double prevAnsInt = 0.;
  double nextAnsInt = 0.;
  double epsabs = 0.;
  double epsrel = 1e-8;
  double abserr = 0.;
  double result = 0.;
  size_t limit = 100;
  bool skipP2 = false;
  struct zeta_params zp;
  //double n000_term = 2. *p.gam *gsl_pow_int(M_SQRTPI,3) *n000_sum( p.q2);
  double n000_term = p.gam *M_PI *n000_sum( p.q2);

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);
  F.function = &integral_zeta_00;
  F.params = &zp;

  while ( relerr_check( prevAnsSum, nextAnsSum)
      || relerr_check( prevAnsInt, nextAnsInt) || (i % 4 != 0) || skipP2 ) {
    //std::cout << "iteration " <<i <<std::endl;
    prevAnsSum = nextAnsSum; prevAnsInt = nextAnsInt; skipP2 = false;
    auto vecCombos = all_combos( i); // get a list of all vector combos for this choice
    if ( vecCombos.size() == 0) { skipP2 = true; }
    for ( auto vecc = vecCombos.begin(); vecc != vecCombos.end(); vecc++ ) {
      auto vecPerms = all_permutations( *vecc);
      for ( auto vecp = vecPerms.begin(); vecp != vecPerms.end(); vecp++ ) {
        //std::cout << "sum vector: "
        //  <<(*vecp)[0] <<", " <<(*vecp)[1] <<", " <<(*vecp)[2];
        zp = full_to_zeta_params( (*vecp)[0], (*vecp)[1], (*vecp)[2], p);
        gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
        nextAnsSum += zp.leadsum;
        nextAnsInt += zp.iphase *result;
        //std::cout << "  |  sum : " << zp.leadsum;
        //std::cout << "  |  int : " << zp.iphase *result << std::endl;
      }
    }
    //std::cout << "skipP2 : " << skipP2 << std::endl;
    //std::cout << "sum: " << prevAnsSum << ", " << nextAnsSum << std::endl;
    //std::cout << "int: " << prevAnsInt << ", " << nextAnsInt << std::endl;
    //std::cout << "sum: " << RELERR << " < "
    //  << abs(prevAnsSum - nextAnsSum)/abs(nextAnsSum) << std::endl;
    //std::cout << "int: " << RELERR << " < "
    //  << abs(prevAnsInt - nextAnsInt)/abs(nextAnsInt) << std::endl;
    i++;
    assert(i < 100);
  }

  //std::cout << "n000 sum : " << n000_sum( p.q2) << std::endl;
  //std::cout << "n000 term: " << n000_term << std::endl;
  return nextAnsSum + nextAnsInt + n000_term;
}

int main(int argc, char** argv)
{

  gsl_set_error_handler( &gsl_to_c_handler ); // set my own error handler

  struct full_params fparams;
  struct zeta_params zparams;
  fparams.dx = 0;
  fparams.dy = 0;
  fparams.dz = 0;
  fparams.q2 = 0.;
  fparams.gam = 1.;
  int nx = 3;
  int ny = 0;
  int nz = 0;
  zparams = full_to_zeta_params( nx,ny,nz, fparams);

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

  //for (int i=1; i<25; i++){
  //  auto vecout = all_combos( i);
  //  std::cout << "i: " <<i << std::endl;
  //  for (auto pvec = vecout.begin(); pvec != vecout.end(); pvec++) {
  //    std::cout << "-- output vector: "
  //    <<(*pvec)[0] <<", " <<(*pvec)[1] <<", " <<(*pvec)[2] << std::endl;
  //  }
  //}

  std::cout << "zeta 00 : " << full_zeta_00 (fparams) << std::endl;
  //std::cout << full_zeta_slow_00 (fparams) << std::endl;
  //std::cout << "dawson: " << erfi( fparams.q2) << std::endl;
  std::cout << "zeta test successful" << std::endl;

  gsl_set_error_handler( NULL );
  return 0;
}
