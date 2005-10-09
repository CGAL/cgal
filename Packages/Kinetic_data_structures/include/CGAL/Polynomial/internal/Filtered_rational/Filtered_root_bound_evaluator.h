#ifndef CGAL_POLYNOMIAL_FILTERED_ROOT_BOUND_EVALUATOR_H
#define CGAL_POLYNOMIAL_FILTERED_ROOT_BOUND_EVALUATOR_H

#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE


template<class Kernel, class M_t = CGAL::Field_tag>
class Filtered_root_bound_evaluator
{
public:
  Filtered_root_bound_evaluator(bool ,
				const Kernel k): rb_(k.interval_traits_object().root_bound_object())  {}
  typedef double result_type;
  typedef typename Kernel::Function argument_type;
  
  result_type operator()(const argument_type& p) const
  {
    Interval_arithmetic_guard iag;
    return rb_(p.interval_function()).sup();
  }
protected:
  typename Kernel::Interval_traits::Root_bound rb_;
};



CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif // CGAL_POLYNOMIAL_ROOT_BOUND_EVALUATOR_H
