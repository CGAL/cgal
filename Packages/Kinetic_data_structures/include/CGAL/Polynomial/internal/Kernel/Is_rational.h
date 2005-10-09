#ifndef CGAL_POLYNOMIAL_INTERNAL_IS_RATIONAL_H
#define CGAL_POLYNOMIAL_INTERNAL_IS_RATIONAL_H

#include <CGAL/Polynomial/basic.h>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

//! Compute the sign after a root.
/*!
  This has specializations for Explicit_roots. 
*/
template <class K>
class Is_rational {
public:
  Is_rational(){  }

  typedef bool result_type;
  typedef typename K::Root argument_type;

  template <class T>
  result_type operator()(const T &v) const {
    return v.is_rational();
  }

  bool operator()(double) const {
    return true;
  }
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
