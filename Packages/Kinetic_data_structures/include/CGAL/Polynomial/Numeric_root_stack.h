#ifndef CGAL_POLYNOMIAL_NUMERIC_SOLVER_H
#define CGAL_POLYNOMIAL_NUMERIC_SOLVER_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Numeric_root_stack_core.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

template <class Solver_traits>
class Numeric_root_stack: public internal::Numeric_root_stack_core<Solver_traits, false> {
  typedef internal::Numeric_root_stack_core<Solver_traits, false> Parent;
public:
  typedef typename Parent::Root Root;
  typedef typename Solver_traits::Function Function;
  Numeric_root_stack(const typename Solver_traits::Function &f, 
			  Root lb, Root ub, 
			  const Solver_traits&k): Parent(f, lb, ub, k){
  }
  Numeric_root_stack(){}
};

CGAL_POLYNOMIAL_END_NAMESPACE
#endif




