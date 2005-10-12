#ifndef CGAL_POLYNOMIAL_INTERNAL_BELOW_AT_ROOT_H
#define CGAL_POLYNOMIAL_INTERNAL_BELOW_AT_ROOT_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Rational/Sign_below_rational.h>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

//! Compute the sign after a root.
/*!
  This has specializations for Explicit_roots. 
*/
template <class R, class K>
class Sign_below: public Sign_below_rational<K>{
  typedef Sign_below_rational<K> P;
public:
  Sign_below(const typename K::Function &p, K k): P(p,k) {
  }
  Sign_below(){}
  using P::operator();
  typename P::result_type operator()(const typename K::Root &v) const {
    CGAL_Polynomial_precondition(0);
    return ZERO;
  }

};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
