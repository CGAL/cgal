#ifndef CGAL_POLYNOMIAL_INTERNAL_EVEN_MULTIPLICITY_H
#define CGAL_POLYNOMIAL_INTERNAL_EVEN_MULTIPLICITY_H

#include <CGAL/Polynomial/basic.h>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

//! Compute the sign after a root.
/*!
  This has specializations for Explicit_roots. 
*/
template <class K>
class Is_even_multiplicity {
public:
  Is_even_multiplicity(){  }

  typedef bool result_type;
  typedef typename K::Root argument_type;

  template <class T>
  result_type operator()(const T &v) const {
    return v.is_even_multiplicity();
  }

  bool operator()(double) const {
    return false;
  }
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
