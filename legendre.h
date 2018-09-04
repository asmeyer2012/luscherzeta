#ifndef __LEGENDRE_H__
#define __LEGENDRE_H__

#include <assert.h>
#include <math.h>
#include <vector>

#include <gsl/gsl_sf_gamma.h>

class poly_term
{
  public:
   poly_term(int zexp, int qexp, double fac);
   poly_term copy() { return poly_term(this->zexp, this->qexp, this->fac); }
   poly_term copy_z_plus_1() { return poly_term(this->zexp+1, this->qexp, this->fac); }
   poly_term copy_q_plus_1() { return poly_term(this->zexp, this->qexp+1, this->fac); }
   std::string as_string();
   void scale(double fac) { this->fac *= fac; }
   void add  (double fac) { this->fac += fac; }
   double evaluate(double z);
   int get_z()   { return this->zexp; }
   int get_q()   { return this->qexp; }
   int get_fac() { return this->fac;  }

  private:
   double fac; // prefactor
   int zexp; // m for z^m
   int qexp; // n for \sqrt{1-z^2}^n
};

poly_term::poly_term(int zexp, int qexp, double fac) {
 this->zexp = zexp;
 this->qexp = qexp;
 this->fac  = fac;
}

std::string poly_term::as_string() {
 stringstream ss;
 if (this->zexp == 0 && this->qexp == 0) { ss << this->fac << "\n"; return ss.str(); }
 else {
   // build qexp term first
   if (this->qexp > 0) {
     if (this->qexp % 2 == 0) { ss << "(1-z^2)^" << this->qexp/2 << " * "; }
     else                     { ss << "(1-z^2)^" << this->qexp   << "/2 * "; }
   }
   if (this->zexp > 0) { ss << "z^"       << this->zexp << " * "; }
   ss << this->fac << "\n";
   return ss.str();
 }
 assert(0);
}

double poly_term::evaluate(double z) {
 return this->fac *pow( sqrt(1.-z*z), this->qexp) *pow(z, this->zexp);
}

// a single legendre polynomial instance, integer order
class associated_legendre_poly
{
  public:
   associated_legendre_poly(std::vector< poly_term> poly= {});
   associated_legendre_poly copy() { return associated_legendre_poly( this->poly); }
   associated_legendre_poly copy_z_plus_1();
   associated_legendre_poly copy_q_plus_1();
   std::string as_string();
   void add(poly_term term);
   void add(associated_legendre_poly alp);
   void scale(double fac);
   double evaluate(double z);
  private:
   std::vector< poly_term> poly;
};

associated_legendre_poly::associated_legendre_poly( std::vector< poly_term> poly) {
  this->poly = poly;
}

associated_legendre_poly associated_legendre_poly::copy_z_plus_1() {
  associated_legendre_poly p0((std::vector< poly_term>) {});
  for (auto term0 = this->poly.begin(); term0 != this->poly.end(); term0++) {
    p0.add( term0->copy_z_plus_1());
  }
  return p0;
}

associated_legendre_poly associated_legendre_poly::copy_q_plus_1() {
  associated_legendre_poly p0((std::vector< poly_term>) {});
  for (auto term0 = this->poly.begin(); term0 != this->poly.end(); term0++) {
    p0.add( term0->copy_q_plus_1());
  }
  return p0;
}

std::string associated_legendre_poly::as_string() {
 stringstream ss;
 int i = 0;
 for (auto term0 = this->poly.begin(); term0 != this->poly.end(); term0++) {
   if (i == 0) { ss << "   "; }
   else        { ss << " + "; }
   ss << term0->as_string();
   i++;
 }
 return ss.str();
}

void associated_legendre_poly::add(poly_term term)
{
  // check for existing terms of the same variety; add to them if they exist
  for (auto term0 = this->poly.begin(); term0 != this->poly.end(); term0++) {
    if (term0->get_z() == term.get_z() && term0->get_q() == term.get_q()) {
     term0->add( term.get_fac());
     return;
    }
  }
  this->poly.push_back( term); // add the unique term
}

void associated_legendre_poly::add(associated_legendre_poly alp)
{
  // add all terms from alp to this one
  for (auto term0 = alp.poly.begin(); term0 != alp.poly.end(); term0++) {
    this->add( *term0);
  }
}

void associated_legendre_poly::scale(double fac)
{
  // scale all terms in polynomial
  for (auto term0 = this->poly.begin(); term0 != this->poly.end(); term0++) {
    term0->scale( fac);
  }
}

double associated_legendre_poly::evaluate(double z)  { 
  double val = 0.;
  for (auto term0 = this->poly.begin(); term0 != this->poly.end(); term0++) {
    val += term0->evaluate(z);
  }
  return val;
}

struct associated_legendre_pair { int l; int m; associated_legendre_poly ascleg; };

// template to build associated legendre polynomials for integer order
class associated_legendre
{
  public:
  associated_legendre(int maxOrder);
  associated_legendre_poly get(int l, int m);
  private:
  std::vector< associated_legendre_pair> poly;
  associated_legendre_poly get_poly(int l, int m);
  int maxOrder;
  void sumrule_lplus1(int l);
  void sumrule_mplus1(int l,int m);
  void sumrule_minusm(int l,int m);
  void increment_max_order();
};

associated_legendre::associated_legendre(int maxOrder = 0) {
  std::cout << "initializing associated_legendre with maxOrder = " << maxOrder << std::endl;
  // fix the 0-order legendre polynomial
  this->maxOrder = 0;
  associated_legendre_poly p0({ poly_term(0,0,1.) });
  associated_legendre_pair alpair = { 0,0,p0 };
  this->poly.push_back(alpair);
  while(maxOrder > this->maxOrder) { this->increment_max_order(); }
  //this->maxOrder = maxOrder;
}

// safe version
associated_legendre_poly associated_legendre::get(int l, int m) {
  assert(abs(m) < l+1);
  while(l > this->maxOrder) { this->increment_max_order(); }
  return this->poly[l*l + l - m].ascleg.copy();
}

// handle version, for manipulation
associated_legendre_poly associated_legendre::get_poly(int l, int m) {
  return this->poly[l*l + l - m].ascleg; // specific order
}

// use the sum rule P_{l+1}^{l+1} = -(2l+1) \sqrt{1-x^2} P_{l}^{l}
void associated_legendre::sumrule_lplus1(int l) {
  associated_legendre_poly pl1 = this->get_poly(l,l).copy_q_plus_1();
  pl1.scale( -(2*l+1.));
  this->poly.push_back( {l+1, l+1, pl1});
}

// use the sum rule (l-m+1) P_{l+1}^{m} = (2l+1)x P_{l}^{m} -(l+m) P_{l-1}^{m}
void associated_legendre::sumrule_mplus1(int l, int m) {
  associated_legendre_poly pl0 = this->get_poly(l,m).copy_z_plus_1();
  associated_legendre_poly pl1;
  if (m < l) { pl1 = this->get_poly(l-1,m).copy(); }
  else       { pl1 = associated_legendre_poly(); } // 0. in poly form
  pl0.scale( (2*l+1)/(l-m+1.));
  pl1.scale(  -(l+m)/(l-m+1.));
  pl0.add( pl1);
  this->poly.push_back( {l+1, m, pl0});
}

// use the sum rule P_{l}^{-m} = Gamma(l-m+1)/Gamma(l+m+1) P_{l}^{m}
void associated_legendre::sumrule_minusm(int l, int m) {
  associated_legendre_poly pl0 = this->get_poly(l+1,-m).copy(); // still old max value of l
  pl0.scale( gsl_sf_gamma( l+m+2) /gsl_sf_gamma( l-m+2));
  this->poly.push_back( {l+1, m, pl0});
}

void associated_legendre::increment_max_order() {
  std::cout << "building legendre polynomials of order " << this->maxOrder+1 << std::endl;
  int l = this->maxOrder;
  for (int m=l+1; m > -l-2; m--) {
    std::cout << "l,m = " << l+1 <<","<< m << std::endl;
    if      (m == l+1) { this->sumrule_lplus1(l);   }
    else if (m >= 0  ) { this->sumrule_mplus1(l,m); }
    else               { this->sumrule_minusm(l,m); }
    //std::cout << this->get_poly(l+1,m).as_string() << std::endl;
  }
  this->maxOrder++;
}

#endif
