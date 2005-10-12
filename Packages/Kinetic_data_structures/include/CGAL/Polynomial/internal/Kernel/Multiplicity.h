#ifndef CGAL_POLYNOMIAL_INTERNAL_MULTIPLICITY_H
#define CGAL_POLYNOMIAL_INTERNAL_MULTIPLICITY_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Rational/Rational_multiplicity.h>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

//! Compute the sign after a root.
/*!
  This has specializations for Explicit_roots. 
*/
template <class K>
class Multiplicity: public Rational_multiplicity<K>{
  typedef Rational_multiplicity<K> P;
public:
  Multiplicity(const typename K::Function &p, K k=K()): P(p,k){  }
  Multiplicity(){}

  using P::operator();
  typename P::result_type operator()(const typename K::Root &v) const {
    CGAL_Polynomial_precondition(0);
    return 1;
  }
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
