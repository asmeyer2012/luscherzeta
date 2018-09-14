#ifndef __GSL_EXCEPTION_H__
#define __GSL_EXCEPTION_H__

#include <exception>
#include <gsl/gsl_errno.h>

using namespace std;

struct gsl_underflow_exception : public exception {
   const char * what () const throw () { return "GSL_UNDERFLOW exception"; }
};

struct gsl_other_exception : public exception {
   const char * what () const throw () { return "GSL_* exception"; }
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

#endif
