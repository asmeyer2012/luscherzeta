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
#include "vector_sorting.h"

#define RELERR 1e-8 // relative error to compute zeta functions to

// complex error function
double erfi( double q2)
{ return 2. *gsl_sf_exp( gsl_pow_int( q2, 2)) *gsl_sf_dawson( q2) /M_SQRTPI; }

// check relative error of two values
bool relerr_check (double prevVal, double nextVal)
{ return ( RELERR *abs( nextVal) < abs( nextVal -prevVal) ); }
bool relerr_check (gsl_complex prevVal, gsl_complex nextVal)
{ return ( RELERR *gsl_complex_abs2( nextVal)
   < gsl_complex_abs2( gsl_complex_sub(nextVal,prevVal)) ); }

// full parameters to get all prefactors correct
struct full_params { int dx; int dy; int dz; double q2; double gam;
 int l; int m; spherical_harmonic * sharm; };

struct zeta_lm_params {
 double r2q2;         // r2 - q2
 double r2;           // r2
 double q2;           // q2 from input
 double gam;          // gam from input
 double ngam2;        // exponent of integral
 gsl_complex leadsum; // leading sum term
 gsl_complex iphase;  // prefactor and phase on integral
 gsl_complex rlYlm;   // factor of r^l Y_{lm}(\hat{r})
 int l;               // total angular momentum quantum number
 int m;               // z-component angular momentum quantum number
 };

struct zeta_lm_params full_to_zeta_lm_params( int inx, int iny, int inz,
  struct full_params in_params )
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
  // included in iphase
  gsl_complex gnlYlm;  // factor of k^l Y_{lm}(\hat{k}) for \vec{k} = 2\pi \hat{\gamma}.\vec{n}

  assert( in_params.l >= abs(in_params.m) ); // also triggers if l < 0

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
    gnx *= 2.*M_PI; gny *= 2.*M_PI; gnz *= 2.*M_PI; // scale by 2\pi for convenience
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
    out_params.rlYlm = gsl_complex_rect( 1./(2.*M_SQRTPI), 0.);
    gnlYlm = out_params.rlYlm;
  } else {
    // Rummukainen-Gottlieb have an extra factor of i^l, should be cancelled by derivative
    out_params.rlYlm = in_params.sharm->evaluate( out_params.l, out_params.m, rx,ry,rz);
    gnlYlm           = in_params.sharm->evaluate( out_params.l, out_params.m, gnx,gny,gnz);
  }

  double r2q2 = r2 - in_params.q2; // argument for leading sum term
  // scale \vec{n} parallel to \vec{d} by \gamma, then square
  // include a factor of \pi^2 for convenience
  out_params.ngam2 = M_PI*M_PI*((gam*gam -1.) *npar2 +n2);
  // == \gamma e^{-i\pi n.d} / 2\sqrt{\pi}, phase for integral
  // nd is always integer, so real
  out_params.iphase = gsl_complex_mul_real( gnlYlm, gam *gsl_pow_int(-1.,int(nd)));
  // other useful terms
  out_params.r2 = r2;
  out_params.r2q2 = r2q2;
  out_params.leadsum = gsl_complex_mul_real( out_params.rlYlm, 1./r2q2 );

  return out_params;
}

// version that subtracts out the divergence at x=0
// only needed for l == 0 since Ylm = 0 at r == 0 for l > 0
double integral_zeta_lm_sub (double x, void * p)
{
  struct zeta_lm_params * params = (struct zeta_lm_params *)p;
  double q2 = (params->q2);
  double ngam2 = (params->ngam2);
  try {
    // deal with roundoff error. thanks Taku
    double y = q2*x -ngam2/x;
    if(fabs(y) > 1e-4) {
      return pow( M_PI/x, 1.5) *(gsl_sf_exp(y) - 1.0);
    }
    else {
      return pow( M_PI/x, 1.5) *y*(1.+y/2*(1.+y/3.*(1.+y/4.*(1.+y/5.*(1.+y/6)))));
    }
    //return pow( M_PI/x, 1.5) *(gsl_sf_exp(q2*x -ngam2/x) - 1.);
  }
  catch ( gsl_underflow_exception& e) {
    return 0.;
  }
}

double integral_zeta_lm(double x, void * p)
{
  struct zeta_lm_params * params = (struct zeta_lm_params *)p;
  double q2 = (params->q2);
  double ngam2 = (params->ngam2);
  int l = (params->l);
  try {
    // (\pi/t)^{3/2} (1/2t)^{l} exp{t q.q - n.n/t}
    return pow( M_PI/x, 1.5) *gsl_pow_int( .5/x, l) *gsl_sf_exp(q2*x -ngam2/x);
  }
  catch ( gsl_underflow_exception& e) {
    return 0.;
  }
}

gsl_complex full_zeta_lm (struct full_params p)
{
  int i = 1;
  gsl_complex prevAns = gsl_complex_rect( 0., 0.);
  gsl_complex nextAns = gsl_complex_rect( 0., 0.);
  double epsabs = 0.;
  double epsrel = 1e-8;
  double abserr = 0.;
  double result = 0.;
  size_t limit = 1000;
  bool skipP2 = false;
  struct zeta_lm_params zp;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);
  F.params = &zp;

  // zero term needs to be handled explicitly
  zp = full_to_zeta_lm_params( 0,0,0, p);
  // Ylm(\vec{r}) / (r2-q2)
  nextAns = gsl_complex_add( nextAns, zp.leadsum);
  // (2\pi)^3 \int_1^\infty dt (\gamma / (4\pi t)^{3/2}) (1/(4\pi)^{1/2}) \delta_{l0} \delta{m0}
  // == \gamma \pi \delta_{l0} \delta{m0}
  if (zp.l == 0) { nextAns = gsl_complex_sub( nextAns, gsl_complex_rect( zp.gam *M_PI, 0.)); }
  // if l > 0, then this term is zero
  // (2\pi)^3 \int_0^1 dt (\gamma / (4\pi t)^{3/2}) (1/(4\pi)^{1/2}) (e^{t q2} - 1)
  // == (\gamma / \sqrt{4\pi}) \int_0^1 dt (\pi /t)^{3/2} (e^{t q2} - 1)
  F.function = &integral_zeta_lm_sub;
  if (zp.l == 0) {
    gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
    nextAns = gsl_complex_add( nextAns, gsl_complex_mul_real( zp.iphase, result));
  }
  // (2\pi)^3 (1/(2\pi)^{3}) Ylm(\vec{r} \int_0^1 dt e^{-t (r2-q2)}
  // == Ylm(\vec{r}) (e^{-(r2-q2)} - 1) / (r2-q2)
  nextAns = gsl_complex_sub( nextAns,
   gsl_complex_mul_real( zp.rlYlm, (1. -gsl_sf_exp(-zp.r2q2)) /zp.r2q2 ));
  // remaining term of \int_1^\infty is an exact cancellation by definition

  prevAns = nextAns;

  F.function = &integral_zeta_lm; // use unsubtracted version for everything else
  while ( i < p.q2 || relerr_check( prevAns, nextAns) || (i % 4 != 0) || skipP2 ) {
    prevAns = nextAns; skipP2 = false;
    auto vecCombos = all_combos( i); // get a list of all vector combos for this choice

    if ( vecCombos.size() == 0) { skipP2 = true; }
    for ( auto vecc = vecCombos.begin(); vecc != vecCombos.end(); vecc++ ) {
      auto vecPerms = all_permutations( *vecc);
      for ( auto vecp = vecPerms.begin(); vecp != vecPerms.end(); vecp++ ) {

        zp = full_to_zeta_lm_params( (*vecp)[0], (*vecp)[1], (*vecp)[2], p);

        if ( zp.r2 < zp.q2) {
          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
          nextAns = gsl_complex_sub( nextAns,
           gsl_complex_mul_real( zp.rlYlm, (1. -gsl_sf_exp(-zp.r2q2)) /zp.r2q2 ));
          nextAns = gsl_complex_add( nextAns, zp.leadsum);
          nextAns = gsl_complex_add( nextAns, gsl_complex_mul_real( zp.iphase, result));
        }
        else {
          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
          nextAns = gsl_complex_add( nextAns,
           gsl_complex_mul_real( zp.rlYlm, gsl_sf_exp(-zp.r2q2) /zp.r2q2 ));
          nextAns = gsl_complex_add( nextAns, gsl_complex_mul_real( zp.iphase, result));
        }
      }
    }

    i++;
    assert(i < 100);
  }

  return nextAns;
}

#endif
