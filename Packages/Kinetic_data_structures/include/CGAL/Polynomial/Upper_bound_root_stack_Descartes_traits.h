#ifndef CGAL_POLYNOMIAL_UPPER_BOUND_ENUMERATOR_DESCARTES_TRAITS_H
#define CGAL_POLYNOMIAL_UPPER_BOUND_ENUMERATOR_DESCARTES_TRAITS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Root_stack_traits_base.h>
#include <CGAL/Polynomial/internal/Rational/Descartes_root_counter.h>
#include <CGAL/Polynomial/internal/Rational/Sturm_root_counter.h>
#include <CGAL/Polynomial/internal/Simple_interval_root.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

template <class Poly>
class Upper_bound_root_stack_Descartes_traits: public internal::Root_stack_traits_base<Poly> {
  typedef internal::Root_stack_traits_base<Poly>  P;
  typedef Upper_bound_root_stack_Descartes_traits<Poly> This;
public:
  typedef POLYNOMIAL_NS::internal::Simple_interval_root<This> Root;

  typedef internal::Descartes_root_counter<This> Root_count;
  Root_count root_count_object(const typename P::Function &f) const {
    return Root_count(f, *this);
  }
  typedef internal::Sturm_root_counter<This> Sturm_root_count;
  Sturm_root_count Sturm_root_count_object(const typename P::Function &f) const {
    return Sturm_root_count(f, *this);
  }
};

CGAL_POLYNOMIAL_END_NAMESPACE

#endif
