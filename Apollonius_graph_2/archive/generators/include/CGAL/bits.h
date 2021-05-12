#ifndef CGAL_BITS_H
#define CGAL_BITS_H

#include <CGAL/basic.h>
#include <cmath>

namespace CGAL {

double log2(double x)
{
  return log10(x) / log10(2.0);
}

unsigned int bits(double x)
{
  CGAL_precondition( static_cast<int>(x) == x );
  if ( x == 0 ) { return 1; }
  return static_cast<unsigned int>(log2( CGAL::abs(x) )) + 1;
}

} //namespace CGAL

#endif // CGAL_BITS_H
