#ifndef CGAL_POLYNOMIAL_INTERNAL_UPPER_BOUND_ENUMERATOR_TRAITS_BASE_H
#define CGAL_POLYNOMIAL_INTERNAL_UPPER_BOUND_ENUMERATOR_TRAITS_BASE_H

#include <CGAL/Polynomial/basic.h>

#include <CGAL/Polynomial/internal/Isolating_interval.h>
#include <CGAL/Polynomial/internal/Rational/Rational_traits_base.h>



CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class Poly>
class Root_stack_traits_base: public Rational_traits_base<Poly> {
private:
  typedef Root_stack_traits_base<Poly> This;
  typedef Rational_traits_base<Poly> P;

public:
  typedef CGAL_POLYNOMIAL_NS::internal::Isolating_interval<typename P::NT> Isolating_interval;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
