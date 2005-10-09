#ifndef CGAL_POLYNOMIAL_INTERNAL_TO_RATIONAL_H
#define CGAL_POLYNOMIAL_INTERNAL_TO_RATIONAL_H

#include <CGAL/Polynomial/basic.h>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

//! Compute the sign after a root.
/*!
  This has specializations for Explicit_roots. 
*/
template <class K>
class To_rational {
public:
  To_rational(){  }

  typedef typename K::NT result_type;
  typedef typename K::Root argument_type;

  template <class T>
  result_type operator()(const T &v) const {
    return v.to_rational();
  }

  double operator()(double v) const {
    return v;
  }
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
