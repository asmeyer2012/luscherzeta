#ifndef __ZETA_H__
#define __ZETA_H__

#include <math.h>
//#include <gsl/gsl_complex.h>      // definitions for complex numbers
//#include <gsl/gsl_complex_math.h> // math with complex numbers
#include <gsl/gsl_integration.h> // GSL integration
#include <gsl/gsl_math.h>        // basic GSL math
#include <gsl/gsl_sf_dawson.h>   // D(x) = \sqrt{\pi}{2} e^{-x^2} Erfi(x)
#include <gsl/gsl_sf_exp.h>      // GSL exponential function
//#include <gsl/gsl_sf_result.h> // for better error handling on GSL, not needed

#define RELERR 1e-8 // relative error to compute zeta functions to

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
  double gam = in_params.gam;
  out_params.q2 = in_params.q2;
  double nd,npar2,r2;

  if (in_params.dx || in_params.dy || in_params.dz) {
    double dx = double(in_params.dx);
    double dy = double(in_params.dy);
    double dz = double(in_params.dz);
    // vector dot products
    double d2 = (dx*dx +dy*dy +dz*dz);
    nd = (nx*dx +ny*dy +nz*dz);
    npar2 = nd*nd/d2;
    // (n_par + d/2)^2/gam^2 +n_perp^2
    r2 = (npar2 +nd +.25*d2)/(gam*gam) +n2 - npar2;

  }
  else {
    // \vec{d} == 0, \gamma == 1.
    nd = 0.;
    npar2 = 0.; // prevent division by 0.
    r2 = n2;
  }

  double r2q2 = r2 - in_params.q2; // argument for leading sum term
  // scale \vec{n} parallel to \vec{d} by \gamma, then square
  // include a factor of \pi^2 for convenience
  out_params.ngam2 = M_PI*M_PI*((gam*gam -1.) *npar2 +n2);
  // == \gamma e^{-i\pi n.d} / 2\sqrt{\pi}, phase for integral
  // nd is always integer, so real
  out_params.iphase = gam /(2.*M_SQRTPI) *gsl_pow_int(-1.,int(nd));
  // other useful terms
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
  nextAns += 1./zp.r2q2/(2.*M_SQRTPI);
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
        nextAns += 1./zp.r2q2/(2.*M_SQRTPI);
        std::cout << "  |  contribution : " << 1./zp.r2q2/(2.*M_SQRTPI) << std::endl;
      }
    }
    std::cout << "skipP2 : " << skipP2 << std::endl;
    std::cout << "prevAns: " << prevAns << std::endl;
    std::cout << "nextAns: " << nextAns << std::endl;
    std::cout << RELERR << " < " << abs(prevAns - nextAns)/abs(nextAns) << std::endl;
    i++;
    assert(i < 1000);
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
  double n000_term = p.gam *M_PI *n000_sum( p.q2);

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);
  F.function = &integral_zeta_00;
  F.params = &zp;

  // add zero term from first sum
  zp = full_to_zeta_params( 0,0,0, p);
  prevAnsSum += zp.leadsum;
  nextAnsSum += zp.leadsum;

  while ( relerr_check( prevAnsSum, nextAnsSum)
      || relerr_check( prevAnsInt, nextAnsInt) || (i % 4 != 0) || skipP2 ) {
    //std::cout << "iteration " <<i <<std::endl;
    prevAnsSum = nextAnsSum; prevAnsInt = nextAnsInt; skipP2 = false;
    auto vecCombos = all_combos( i); // get a list of all vector combos for this choice
    if ( vecCombos.size() == 0) { skipP2 = true; }
    for ( auto vecc = vecCombos.begin(); vecc != vecCombos.end(); vecc++ ) {
      auto vecPerms = all_permutations( *vecc);
      for ( auto vecp = vecPerms.begin(); vecp != vecPerms.end(); vecp++ ) {
        zp = full_to_zeta_params( (*vecp)[0], (*vecp)[1], (*vecp)[2], p);
        gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
        nextAnsSum += zp.leadsum;
        nextAnsInt += zp.iphase *result;
        //std::cout << "sum vector: "
        //  <<(*vecp)[0] <<", " <<(*vecp)[1] <<", " <<(*vecp)[2];
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

  //std::cout << "sum term: " << nextAnsSum << std::endl;
  //std::cout << "int term: " << nextAnsInt << std::endl;
  //std::cout << "000 term: " << n000_term << std::endl;
  return nextAnsSum + nextAnsInt + n000_term;
}

#endif