#ifndef __ZETA_H__
#define __ZETA_H__

#include <assert.h>
#include <math.h>
#include <gsl/gsl_complex.h>      // definitions for complex numbers
#include <gsl/gsl_complex_math.h> // math with complex numbers
#include <gsl/gsl_integration.h> // GSL integration
#include <gsl/gsl_math.h>        // basic GSL math
#include <gsl/gsl_sf_dawson.h>   // D(x) = \sqrt{\pi}{2} e^{-x^2} Erfi(x)
#include <gsl/gsl_sf_exp.h>      // GSL exponential function
//#include <gsl/gsl_sf_result.h> // for better error handling on GSL, not needed

#include "spherical.h"

#define RELERR 1e-8 // relative error to compute zeta functions to

// complex error function
double erfi( double q2)
{ return 2. *gsl_sf_exp( gsl_pow_int( q2, 2)) *gsl_sf_dawson( q2) /M_SQRTPI; }

// check relative error of two values
bool relerr_check (double prevVal, double nextVal)
{ return ( RELERR *abs( nextVal) < abs( nextVal -prevVal) ); }

// full parameters to get all prefactors correct
struct full_params { int dx; int dy; int dz; double q2; double gam; int l; int m; };
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

// integration method, faster and convergent
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
      }
    }
    i++;
    assert(i < 100);
  }

  return nextAnsSum + nextAnsInt + n000_term;
}

// just compute the parameters that are needed
struct zeta_med_params {
 double r2q2;        // r2 - q2
 double r2;          // r2
 double q2;          // q2 from input
 double gam;         // gam from input
 double ngam2;       // exponent of integral
 double leadsum;     // leading sum term
 double iphase;      // prefactor and phase on integral
 };

// \sum_l (q^2)^l /((2l-1) \Gamma(l+1) ) == (\sqrt{q^2} Erfi{\sqrt{q^2}} - e^{q^2})
struct zeta_med_params full_to_zeta_med_params( int inx, int iny, int inz,
  const struct full_params in_params )
{
  struct zeta_med_params out_params;
  double nx = double(inx);
  double ny = double(iny);
  double nz = double(inz);
  double n2 = (nx*nx +ny*ny +nz*nz);
  double gam = in_params.gam;
  out_params.q2 = in_params.q2;
  out_params.gam = in_params.gam;
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
  out_params.r2 = r2;
  out_params.r2q2 = r2q2;
  out_params.leadsum = 1. /(2.*M_SQRTPI*r2q2);
  return out_params;
}

// version that subtracts out the divergence at x=0
double integral_zeta_00_sub (double x, void * p)
{
  struct zeta_med_params * params = (struct zeta_med_params *)p;
  double q2 = (params->q2);
  double ngam2 = (params->ngam2);
  // underflows are okay, just give 0. good enough.
  try {
    return pow( M_PI/x, 1.5) *(gsl_sf_exp(q2*x -ngam2/x) - 1.);
  }
  catch ( gsl_underflow_exception& e) {
    return 0.;
  }
}

double integral_zeta_00_med(double x, void * p)
{
  struct zeta_med_params * params = (struct zeta_med_params *)p;
  double q2 = (params->q2);
  double ngam2 = (params->ngam2);
  // underflows are okay, just give 0. good enough.
  try {
    return pow( M_PI/x, 1.5) *gsl_sf_exp(q2*x -ngam2/x); // exp{t q.q - n.n/t}
  }
  catch ( gsl_underflow_exception& e) {
    return 0.;
  }
}

// slower method that should work for lm nonzero
double full_zeta_med_00 (struct full_params p)
{
  int i = 1;
  double prevAns = 0.;
  double nextAns = 0.;
  double epsabs = 0.;
  double epsrel = 1e-6;
  double abserr = 0.;
  double result = 0.;
  size_t limit = 1000;
  bool skipP2 = false;
  struct zeta_med_params zp;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);
  F.params = &zp;

  // zero term needs to be handled explicitly
  zp = full_to_zeta_med_params( 0,0,0, p);
  // (4\pi)^{-1/2} / (r2-q2)
  nextAns += zp.leadsum;
  // (2\pi)^3 \int_1^\infty dt (\gamma / (4\pi t)^{3/2}) (1/(4\pi)^{1/2})
  // == \gamma \pi
  nextAns -= zp.gam *M_PI;
  // (2\pi)^3 \int_0^1 dt (\gamma / (4\pi t)^{3/2}) (1/(4\pi)^{1/2}) (e^{t q2} - 1)
  // == (\gamma / \sqrt{4\pi}) \int_0^1 dt (\pi /t)^{3/2} (e^{t q2} - 1)
  F.function = &integral_zeta_00_sub;
  gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
  nextAns += zp.iphase *result;
  // (2\pi)^3 (1/(2\pi)^{3}) (4\pi)^{-1/2} \int_0^1 dt e^{-t (r2-q2)}
  // == (4\pi)^{-1/2} (e^{-(r2-q2)} - 1) / (r2-q2)
  nextAns -= (1. -gsl_sf_exp(-zp.r2q2)) /zp.r2q2 /(2.*M_SQRTPI);
  // remaining term of \int_1^\infty is an exact cancellation by definition

  prevAns = nextAns;

  F.function = &integral_zeta_00_med; // use unsubtracted version for everything else
  while ( i < p.q2 || relerr_check( prevAns, nextAns) || (i % 4 != 0) || skipP2 ) {
    prevAns = nextAns; skipP2 = false;
    auto vecCombos = all_combos( i); // get a list of all vector combos for this choice

    if ( vecCombos.size() == 0) { skipP2 = true; }
    for ( auto vecc = vecCombos.begin(); vecc != vecCombos.end(); vecc++ ) {
      auto vecPerms = all_permutations( *vecc);
      for ( auto vecp = vecPerms.begin(); vecp != vecPerms.end(); vecp++ ) {

        zp = full_to_zeta_med_params( (*vecp)[0], (*vecp)[1], (*vecp)[2], p);

        if ( zp.r2 < p.q2) {
          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
          nextAns -= (1. -gsl_sf_exp(-zp.r2q2)) /zp.r2q2 /(2.*M_SQRTPI);
          nextAns += zp.leadsum;
          nextAns += zp.iphase *result;
        }
        else {
          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
          nextAns += gsl_sf_exp(-zp.r2q2) /zp.r2q2 /(2.*M_SQRTPI);
          nextAns += zp.iphase *result;
        }
      }
    }

    i++;
    assert(i < 100);
  }

  return nextAns;
}

// just compute the parameters that are needed
struct zeta_lm_params {
 double r2q2;        // r2 - q2
 double r2;          // r2
 double q2;          // q2 from input
 double gam;         // gam from input
 double ngam2;       // exponent of integral
 double leadsum;     // leading sum term
 double iphase;      // prefactor and phase on integral
 gsl_complex rlYlm;  // factor of r^l Y_{lm}(\hat{r})
 gsl_complex gnlYlm; // factor of k^l Y_{lm}(\hat{k}) for \vec{k} = 2\pi \hat{\gamma}.\vec{n}
 int l;              // total angular momentum quantum number
 int m;              // z-component angular momentum quantum number
 };

struct zeta_lm_params full_to_zeta_lm_params( int inx, int iny, int inz,
  spherical_harmonic * sharm, const struct full_params in_params )
{
  struct zeta_lm_params out_params;
  double nx = double(inx);
  double ny = double(iny);
  double nz = double(inz);
  double n2 = (nx*nx +ny*ny +nz*nz);
  double gam = in_params.gam;
  out_params.q2 = in_params.q2;
  out_params.gam = in_params.gam;
  out_params.l = in_params.l;
  out_params.m = in_params.m;
  double nd,npar2,r2;
  double gnx,gny,gnz; // \hat{\gamma}.\vec{n}
  double rx,ry,rz;

  assert( in_params.l < abs(in_params.m) ); // also triggers if l < 0

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

    // compute \hat{\gamma}.\vec{n}
    double dfac = (gam-1.)*nd/d2;
    gnx = dfac*dx + nx; gny = dfac*dy + ny; gnz = dfac*dz + nz;
    // compute \vec{r} = \hat{\gamma}^{-1}.(\vec{n} + \vec{d}/2)
    dfac = ((1.-gam)*nd/d2 + .5) /gam;
    rx = dfac*dx + nx; ry = dfac*dy + ny; rz = dfac*dz + nz;
  }
  else {
    // \vec{d} == 0, \gamma == 1.
    nd = 0.;
    npar2 = 0.; // prevent division by 0.
    r2 = n2;
    gnx = nx; gny = ny; gnz = nz;
    rx = nx; ry = ny; rz = nz;
  }

  // angular momentum
  if ( in_params.l == 0 ) {
    out_params.rlYlm = gsl_complex_polar( 1./(2.*M_SQRTPI), 0.);
    out_params.gnlYlm = out_params.rlYlm;
  } else {
    out_params.rlYlm  = sharm->evaluate( out_params.l, out_params.m, rx,ry,rz);
    out_params.gnlYlm = sharm->evaluate( out_params.l, out_params.m, gnx,gny,gnz);
  }

  double r2q2 = r2 - in_params.q2; // argument for leading sum term
  // scale \vec{n} parallel to \vec{d} by \gamma, then square
  // include a factor of \pi^2 for convenience
  out_params.ngam2 = M_PI*M_PI*((gam*gam -1.) *npar2 +n2);
  // == \gamma e^{-i\pi n.d} / 2\sqrt{\pi}, phase for integral
  // nd is always integer, so real
  out_params.iphase = gam /(2.*M_SQRTPI) *gsl_pow_int(-1.,int(nd));
  // other useful terms
  out_params.r2 = r2;
  out_params.r2q2 = r2q2;
  out_params.leadsum = 1. /(2.*M_SQRTPI*r2q2);

  return out_params;
}

//// version that subtracts out the divergence at x=0
//double integral_zeta_00_sub (double x, void * p)
//{
//  struct zeta_med_params * params = (struct zeta_med_params *)p;
//  double q2 = (params->q2);
//  double ngam2 = (params->ngam2);
//  // underflows are okay, just give 0. good enough.
//  try {
//    return pow( M_PI/x, 1.5) *(gsl_sf_exp(q2*x -ngam2/x) - 1.);
//  }
//  catch ( gsl_underflow_exception& e) {
//    return 0.;
//  }
//}
//
//double integral_zeta_00_med(double x, void * p)
//{
//  struct zeta_med_params * params = (struct zeta_med_params *)p;
//  double q2 = (params->q2);
//  double ngam2 = (params->ngam2);
//  // underflows are okay, just give 0. good enough.
//  try {
//    return pow( M_PI/x, 1.5) *gsl_sf_exp(q2*x -ngam2/x); // exp{t q.q - n.n/t}
//  }
//  catch ( gsl_underflow_exception& e) {
//    return 0.;
//  }
//}
//
//// slower method that should work for lm nonzero
//double full_zeta_med_00 (struct full_params p)
//{
//  int i = 1;
//  double prevAns = 0.;
//  double nextAns = 0.;
//  double epsabs = 0.;
//  double epsrel = 1e-6;
//  double abserr = 0.;
//  double result = 0.;
//  size_t limit = 1000;
//  bool skipP2 = false;
//  struct zeta_med_params zp;
//
//  gsl_function F;
//  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);
//  F.params = &zp;
//
//  // zero term needs to be handled explicitly
//  zp = full_to_zeta_med_params( 0,0,0, p);
//  // (4\pi)^{-1/2} / (r2-q2)
//  nextAns += zp.leadsum;
//  // (2\pi)^3 \int_1^\infty dt (\gamma / (4\pi t)^{3/2}) (1/(4\pi)^{1/2})
//  // == \gamma \pi
//  nextAns -= zp.gam *M_PI;
//  // (2\pi)^3 \int_0^1 dt (\gamma / (4\pi t)^{3/2}) (1/(4\pi)^{1/2}) (e^{t q2} - 1)
//  // == (\gamma / \sqrt{4\pi}) \int_0^1 dt (\pi /t)^{3/2} (e^{t q2} - 1)
//  F.function = &integral_zeta_00_sub;
//  gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
//  nextAns += zp.iphase *result;
//  // (2\pi)^3 (1/(2\pi)^{3}) (4\pi)^{-1/2} \int_0^1 dt e^{-t (r2-q2)}
//  // == (4\pi)^{-1/2} (e^{-(r2-q2)} - 1) / (r2-q2)
//  nextAns -= (1. -gsl_sf_exp(-zp.r2q2)) /zp.r2q2 /(2.*M_SQRTPI);
//  // remaining term of \int_1^\infty is an exact cancellation by definition
//
//  prevAns = nextAns;
//
//  F.function = &integral_zeta_00_med; // use unsubtracted version for everything else
//  while ( i < p.q2 || relerr_check( prevAns, nextAns) || (i % 4 != 0) || skipP2 ) {
//    prevAns = nextAns; skipP2 = false;
//    auto vecCombos = all_combos( i); // get a list of all vector combos for this choice
//
//    if ( vecCombos.size() == 0) { skipP2 = true; }
//    for ( auto vecc = vecCombos.begin(); vecc != vecCombos.end(); vecc++ ) {
//      auto vecPerms = all_permutations( *vecc);
//      for ( auto vecp = vecPerms.begin(); vecp != vecPerms.end(); vecp++ ) {
//
//        zp = full_to_zeta_med_params( (*vecp)[0], (*vecp)[1], (*vecp)[2], p);
//
//        if ( zp.r2 < p.q2) {
//          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
//          nextAns -= (1. -gsl_sf_exp(-zp.r2q2)) /zp.r2q2 /(2.*M_SQRTPI);
//          nextAns += zp.leadsum;
//          nextAns += zp.iphase *result;
//        }
//        else {
//          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
//          nextAns += gsl_sf_exp(-zp.r2q2) /zp.r2q2 /(2.*M_SQRTPI);
//          nextAns += zp.iphase *result;
//        }
//      }
//    }
//
//    i++;
//    assert(i < 100);
//  }
//
//  return nextAns;
//}

#endif
