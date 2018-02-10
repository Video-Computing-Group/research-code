#ifndef __MathFunctions_h
#define __MathFunctions_h

#include <assert.h>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#ifdef _MSC_VER
#include ".\mathplus.h"
#endif

#ifndef exp2
inline double exp2(double x) throw()
{ 
  return exp(x*log(2.0));
}
#endif

#ifndef log2
inline double log2(double x) throw() 
{
  return log(x)/log(2.0);
}
#endif

///
inline double lChoose( double k , double n )
{
  if( k > 0 && k < n )
    return (lgamma( n + 1 ) - lgamma( k + 1 ) - lgamma( n - k + 1 ))/log(2.0);
  else
    return 0;
}

///
inline double l2gamma(double x )
{
  return lgamma(x) / log(2.0);
}

///
inline double digamma(double x)
{
  return ( log(x) - 1/(2*x) - 1/(12*pow(x,2)) +
	   1/(120*pow(x,4)) - 1/(252*pow(x,6)));
}

///
inline double digamma1(double x)
{
  return (1/x + 1/(2*pow(x,2) + 1/(6*pow(x,3)) -
	   1/(30*pow(x,5)) + 1/(42*pow(x,7)) ));
}

///
inline double
GaussPDF(double x, double mu, double sigma2 )
{
  double z = (x -mu)/sqrt(sigma2);
  return 0.5*( 1 + erf(z / M_SQRT2 ));
}

///
inline double
Gausspdf(double x, double mu, double sigma2 )
{
  double z = (x -mu);
  return 1/sqrt(2*M_PI*sigma2)*exp(-0.5*z*z/sigma2);
}

///
/*
inline double
AddLog(double x, double y )
{
  if( x == -HUGE_VAL ) return y;
  if( y == -HUGE_VAL ) return x;
  
  double z = max(x,y);
  return z + log(exp(x-z) + exp(y-z));
}
*/


///
inline double
SubLog2(double x, double y )
{
  if ( x >= y  )
    return x + log2(1 - exp2(y-x));
  else
    return y + log2(exp2(x-y) - 1);
}

///
inline double
SubLog(double x, double y )
{
  if ( x >= y )
    return x + log(1 - exp(y-x));
  else
    return y + log(exp(x-y) - 1);
}

// If x = log2 a, y = log2 b, returns log2 |a-b|
inline double
AbsSubLog2(double x, double y )
{
  if (x == -HUGE_VAL && y == -HUGE_VAL) {
    //log2 |2^-HUGE_VAL - 2^-HUGE_VAL| = log2 |0-0| = log2(0) = -HUGE_VAL
    return -HUGE_VAL;
  }
  else if ( x >= y  )
    return x + log2(1 - exp2(y-x));
  else
    return y + log2(1 - exp2(x-y));
}

// If x = log a, y = log b, returns log |a-b|
inline double
AbsSubLog(double x, double y )
{
  if (x == -HUGE_VAL && y == -HUGE_VAL) {
    //log |e^-HUGE_VAL - e^-HUGE_VAL| = log |0-0| = log(0) = -HUGE_VAL
    return -HUGE_VAL;
  }
  else if ( x >= y )
    return x + log(1.0 - exp(y-x));
  else
    return y + log(1.0 - exp(x-y));
}

inline double
AbsSubLogFactor(double x, double y, double factor)
{
  if (x == -HUGE_VAL && y == -HUGE_VAL) {
    //log |e^-HUGE_VAL - e^-HUGE_VAL| = log |0-0| = log(0) = -HUGE_VAL
    return -HUGE_VAL;
  }
  else if ( (x/factor) >= (y/factor) )
    return x + factor*log(1.0 - exp((y-x)/factor));
  else
    return y + factor*log(1.0 - exp((x-y)/factor));
}

inline double
AddLog(double x, double y )
{
  if( x == -HUGE_VAL ) return y;
  if( y == -HUGE_VAL ) return x;
  
  if (x>=y) {
    return x + log(1.0 + exp(y-x));
  }
  else {
    return y + log(1.0 + exp(x-y));
  }
}

inline double
AddLogFactor(double x, double y, double factor)
{
  if( x == -HUGE_VAL ) return y;
  if( y == -HUGE_VAL ) return x;

  if ((x/factor)>=(y/factor)) { 
    return x + factor * log(1.0 + exp((y-x)/factor));
  }
  else {
    return y + factor * log(1.0 + exp((x-y)/factor));
  }
}

///
inline double
AddLog2(double x, double y )
{
  if( x == -HUGE_VAL ) return y;
  if( y == -HUGE_VAL ) return x;

  double z;
  if( x > y )
    z = y - x;
  else
  {
    z = x - y;
    x = y;
  }
  return x + log2(1 + exp2(z));
}

#define EPS 0.000000001
#define UNDEF_VAL    -1000000000
#define UNDEF_PROB   1.0

#endif
