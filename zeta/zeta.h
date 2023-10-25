#ifndef __ZETA_H__
#define __ZETA_H__

#include <assert.h>
#include <math.h>
#include <limits>
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

// check relative error of two values
bool relerr_check (double prevVal, double nextVal)
{ return ( RELERR *abs( nextVal) < abs( nextVal -prevVal) ); }

bool relerr_check (gsl_complex prevVal, gsl_complex nextVal)
{ return ( RELERR *gsl_complex_abs2( nextVal)
   < gsl_complex_abs2( gsl_complex_sub(nextVal,prevVal)) ); }

// full input parameters to get all prefactors correct
struct full_params {
  // center of mass frame momentum units, shifted for particle mass difference
  double sx;
  double sy;
  double sz;
  // dimensionless momentum and lorentz factor
  double u2;
  double gamma;
  // spherical harmonic total and z-component angular momentum
  int l;
  int m;
  spherical_harmonic * sharm;
  };

// full_params values converted to commonly used quantities
struct zeta_lm_params {
  double z2;           // z2 = ( (n_parallel + s/2)^2 / gamma^2 ) +n_perp^2
  double u2;           // u2 from input
  double z2u2;         // z2 - u2
  double gamma;        // lorentz factor, gamma = EL/ECM = EL/sqrt( EL^2-P^2 )
  double ngam2;        // scaled exponent of integral: pi^2 ( (gamma^2-1) n_parallel^2 +n^2 )
  gsl_complex leadsum; // leading sum term
  gsl_complex iphase;  // prefactor and phase on integral
  gsl_complex rlYlm;   // factor of r^l Y_{lm}(\hat{r})
  int l;               // total angular momentum quantum number
  int m;               // z-component angular momentum quantum number
  };

struct zeta_lm_params full_to_zeta_lm_params(
  int nx0, int ny0, int nz0,
  struct full_params in_params )
{
  struct zeta_lm_params out_params;
  double nx = double(nx0);
  double ny = double(ny0);
  double nz = double(nz0);
  double n2 = (nx*nx +ny*ny +nz*nz);
  double gamma = in_params.gamma;

  assert( in_params.l >= abs(in_params.m) ); // also triggers if l < 0

  // some quick copy assignments
  out_params.u2 = in_params.u2;
  out_params.gamma = in_params.gamma;
  out_params.l = in_params.l;
  out_params.m = in_params.m;

  // intermediates
  //double ns, npar2, z2;
  //double gnx, gny, gnz; // \hat{\gamma}.\vec{n}
  //double zx, zy, zz; // \vec{n} +\gamma^{-1} ( 1/2 + (1-\gamma) \vec{n}\cdot\vec{s} s^{-2} )\vec{s}

  // included in iphase
  // factor of k^l Y_{lm}(\hat{k}) for \vec{k} = 2\pi \hat{\gamma}.\vec{n}
  gsl_complex gnlYlm;

  //if (in_params.dx || in_params.dy || in_params.dz) {
  double sx = in_params.sx;
  double sy = in_params.sy;
  double sz = in_params.sz;

  // vector dot products
  double s2 = (sx*sx +sy*sy +sz*sz);
  double ns = (nx*sx +ny*sy +nz*sz);

  double npar2, z2, sfac0;
  if ( s2 > 1.0e-15 ) {
    // \vec{n} parallel to \vec{s}, magnitude squared
    npar2 = ns*ns /s2;

    // (n_par + s/2)^2 /gam^2 +n_perp^2
    z2 = (npar2 +ns +.25*s2) /(gamma*gamma) +n2 -npar2;

    // this factor shows up a lot in z, w
    // (1-\gamma) \vec{n}\cdot\vec{s} /s^2
    sfac0 = (1. -gamma) *ns /s2;

  } else {
    // prevent division by 0 from s2
    npar2 = sfac0 = 0.;
    z2 = n2;
  }

  // lorentz boosted \vec{n}: \gamma \vec{n}_\parallel +\vec{n}_\perp
  // \vec{n} + (1-\gamma) \vec{n}\cdot\vec{s} s^{-2} \vec{s}
  double gnx = nx -sfac0*sx;
  double gny = ny -sfac0*sy;
  double gnz = nz -sfac0*sz;
  // scale by 2\pi for convenience
  gnx *= 2. *M_PI;
  gny *= 2. *M_PI;
  gnz *= 2. *M_PI;

  // this factor shows up in z
  // \gamma^{-1} ( 1/2 + (1-\gamma) \vec{n}\cdot\vec{s} s^{-2} )
  // \vec{z} = \vec{n} +\gamma^{-1} ( 1/2 + (1-\gamma) \vec{n}\cdot\vec{s} s^{-2} )\vec{s}
  double sfac1 = (.5 +sfac0) /gamma;
  double zx = sfac1 *sx +nx;
  double zy = sfac1 *sy +ny;
  double zz = sfac1 *sz +nz;

  // angular momentum
  if ( in_params.l == 0 ) {
    out_params.rlYlm = gsl_complex_rect( 1./(2.*M_SQRTPI), 0.);
    gnlYlm = out_params.rlYlm;
  } else {
    // Rummukainen-Gottlieb have an extra factor of i^l, should be cancelled by derivative
    out_params.rlYlm = in_params.sharm->evaluate( out_params.l, out_params.m, zx , zy , zz );
    gnlYlm           = in_params.sharm->evaluate( out_params.l, out_params.m, gnx, gny, gnz);
  }

  double z2u2 = z2 - in_params.u2; // argument for leading sum term

  // scale \vec{n} parallel to \vec{d} by \gamma, then square
  // include a factor of \pi^2 for convenience
  out_params.ngam2 = M_PI*M_PI*((gamma*gamma -1.) *npar2 +n2);

  // == \gamma e^{-i\pi n.d} / 2\sqrt{\pi}, phase for integral
  // nd is always integer, so real
  //out_params.iphase = gsl_complex_mul_real( gnlYlm, gam *gsl_pow_int(-1.,int(nd)));
  gsl_complex exp_exponent = gsl_complex_rect( 0, M_PI *ns);
  gsl_complex exp_product  = gsl_complex_mul( gnlYlm, gsl_complex_exp( exp_exponent ));
  out_params.iphase = gsl_complex_mul_real( exp_product, gamma);

  // other useful terms
  out_params.z2 = z2;
  out_params.z2u2 = z2u2;
  out_params.leadsum = gsl_complex_mul_real( out_params.rlYlm, 1./z2u2 );

  return out_params;
}

// version that subtracts out the divergence at x=0
// only needed for l == 0 since Ylm = 0 at r == 0 for l > 0
double integral_zeta_lm_sub (double x, void * p)
{
  struct zeta_lm_params * params = (struct zeta_lm_params *)p;
  double u2 = (params->u2);
  double ngam2 = (params->ngam2);
  try {
    // deal with roundoff error. thanks Taku
    double y = u2*x -ngam2/x;
    if(fabs(y) > 1e-4) {
      return pow( M_PI/x, 1.5) *(gsl_sf_exp(y) - 1.0);
    }
    else {
      return pow( M_PI/x, 1.5) *y*(1.+y/2*(1.+y/3.*(1.+y/4.*(1.+y/5.*(1.+y/6)))));
    }
    //return pow( M_PI/x, 1.5) *(gsl_sf_exp(u2*x -ngam2/x) - 1.);
  }
  catch ( gsl_underflow_exception& e) {
    return 0.;
  }
}

double integral_zeta_lm(double x, void * p)
{
  struct zeta_lm_params * params = (struct zeta_lm_params *)p;
  double u2 = (params->u2);
  double ngam2 = (params->ngam2);
  int l = (params->l);
  try {
    // (\pi/t)^{3/2} (1/2t)^{l} exp{t q.q - n.n/t}
    return pow( M_PI/x, 1.5) *gsl_pow_int( .5/x, l) *gsl_sf_exp( u2*x -ngam2/x);
  }
  catch ( gsl_underflow_exception& e) {
    return 0.;
  }
}

gsl_complex full_zeta_lm (struct full_params p)
{
  int i = 1; // iteration counter
  gsl_complex prevAns = gsl_complex_rect( 0., 0.);
  gsl_complex nextAns = gsl_complex_rect( 0., 0.);
  double epsabs = 0.;
  double epsrel = 1e-8;
  double abserr = 0.;
  double result = 0.;
  size_t limit = 1000;
  struct zeta_lm_params zp;

  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);
  F.params = &zp;

  // (nx,ny,nz) = 0 term needs to be handled explicitly
  // Ylm(\vec{r}) / (z2-u2)
  zp = full_to_zeta_lm_params( 0,0,0, p);
  nextAns = gsl_complex_add( nextAns, zp.leadsum);

  // (2\pi)^3 \int_1^\infty dt (\gamma / (4\pi t)^{3/2}) (1/(4\pi)^{1/2}) \delta_{l0} \delta{m0}
  // == \gamma \pi \delta_{l0} \delta{m0}
  if (zp.l == 0) {
    nextAns = gsl_complex_sub( nextAns, gsl_complex_rect( zp.gamma *M_PI, 0.));
  }

  // if l > 0, then this term is zero
  // (2\pi)^3 \int_0^1 dt (\gamma / (4\pi t)^{3/2}) (1/(4\pi)^{1/2}) (e^{t u2} - 1)
  // == (\gamma / \sqrt{4\pi}) \int_0^1 dt (\pi /t)^{3/2} (e^{t u2} - 1)
  F.function = &integral_zeta_lm_sub;
  if (zp.l == 0) {
    gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
    nextAns = gsl_complex_add( nextAns, gsl_complex_mul_real( zp.iphase, result));
  }

  // (2\pi)^3 (1/(2\pi)^{3}) Ylm(\vec{r} \int_0^1 dt e^{-t (z2-u2)}
  // == Ylm(\vec{r}) (e^{-(z2-u2)} - 1) / (z2-u2)
  nextAns = gsl_complex_sub( nextAns,
   gsl_complex_mul_real( zp.rlYlm, (1. -gsl_sf_exp(-zp.z2u2)) /zp.z2u2 ));

  // remaining term of \int_1^\infty is an exact cancellation by definition

  prevAns = nextAns;

  // sum the contributions of the remaining (nx,ny,nz) != 0 tuples
  // use unsubtracted version for everything else
  F.function = &integral_zeta_lm;
  while ( i < p.u2 || relerr_check( prevAns, nextAns) || (i % 4 != 0) ) {
    prevAns = nextAns;
    auto vecCombos = all_combos( i); // get a list of all vector combos for this choice

    for ( auto vecc = vecCombos.begin(); vecc != vecCombos.end(); vecc++ ) {
      auto vecPerms = all_permutations( *vecc);
      for ( auto vecp = vecPerms.begin(); vecp != vecPerms.end(); vecp++ ) {

        zp = full_to_zeta_lm_params( (*vecp)[0], (*vecp)[1], (*vecp)[2], p);

        // qag integration over: \int_0^1 dt (\pi/t)^{3/2} (1/2t)^{l} exp{t q.q - n.n/t}
        if ( zp.z2 < zp.u2) {
          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
          nextAns = gsl_complex_sub( nextAns,
           gsl_complex_mul_real( zp.rlYlm, (1. -gsl_sf_exp(-zp.z2u2)) /zp.z2u2 ));
          nextAns = gsl_complex_add( nextAns, zp.leadsum);
          nextAns = gsl_complex_add( nextAns, gsl_complex_mul_real( zp.iphase, result));
        }
        else {
          gsl_integration_qag(&F, 0., 1., epsabs, epsrel, limit, 1, w, &result, &abserr);
          nextAns = gsl_complex_add( nextAns,
           gsl_complex_mul_real( zp.rlYlm, gsl_sf_exp(-zp.z2u2) /zp.z2u2 ));
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
