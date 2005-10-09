#ifndef POLYNOMIAL_NUMERIC_SOLVER_CORE_H
#define POLYNOMIAL_NUMERIC_SOLVER_CORE_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>
#include <CGAL/Polynomial/Polynomial.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class Solver_traits, bool CLEAN>
class Numeric_root_stack_core {

  template <class Fn>
  void initialize(const Fn &f, double lb, double ub) {
    std::vector<double> c(f.degree()+1);
    for (unsigned int i=0; i<= f.degree(); ++i){
      c[i]= to_double(f[i]);
    }
    if (CLEAN) {
      polynomial_compute_cleaned_roots(&*c.begin(), &*c.begin()+ f.degree()+1, lb, ub, roots_);
    } else {
      polynomial_compute_roots(&*c.begin(), &*c.begin()+ f.degree()+1, lb, ub, roots_);
    }
  }
  void initialize(const Polynomial<double> &f, double lb, double ub) {
    if (CLEAN) {
      polynomial_compute_cleaned_roots(&*f.begin(), &*f.begin()+ f.degree()+1, lb, ub, roots_);
    } else {
      polynomial_compute_roots(&*f.begin(), &*f.begin()+ f.degree()+1, lb, ub, roots_);
    }
  }
public:
  typedef double Root;
  typedef typename Solver_traits::Function Function;
  typedef Solver_traits Traits;
  Numeric_root_stack_core(const typename Solver_traits::Function &f, Root lb, Root ub, const Solver_traits&){
    initialize(f, lb, ub);
  }
  
  Numeric_root_stack_core(){};

  void pop() {
    Polynomial_precondition(!roots_.empty());
    roots_.pop_back();
  }

  const Root& top() const {
    Polynomial_precondition(!roots_.empty());
    return roots_.back();
  }

  bool empty() const {
    return roots_.empty();
  }

protected:
  std::vector<double> roots_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif




