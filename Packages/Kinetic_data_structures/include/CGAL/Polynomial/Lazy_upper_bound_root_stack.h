#ifndef CGAL_POLYNOMIAL_LAZY_UPPER_BOUND_ROOT_STACK_H
#define CGAL_POLYNOMIAL_LAZY_UPPER_BOUND_ROOT_STACK_H
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/Tools/Ref_counted.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

/*
  Rep stores unisolated part and has a single pointer to the place to put the next isolated root (function stays in rep).

  Roots check themselves first, then isolated, then rep.


*/


template <class Traits_t, class Rep>
class Lazy_upper_bound_root_stack_root {
  typedef Lazy_upper_bound_root_stack_root<Traits_t, Rep> This;

  int compare(const This &o) const {
    if (iv_.inf() > o.iv_.sup()) return 1;
    else if (iv_.sup() < o.iv_.inf()) return -1;
    else if (iv_.is_point
  }
public:
  template <class NT>
  Lazy_upper_bound_root_stack_root(const NT &nt): iv_(to_interval(nt)), root_(new Rep::Shared_root(nt)) {}
  Lazy_upper_bound_root_stack_root(double d): iv_(d) {}
  
  bool operator==(const This &o) const {
    return compare(o)==0;
  }
  bool operator<(const This &o) const {
    return compare(o)==-1;
  }
  bool operator>(const This &o) const {
    return compare(o)==1;
  }
  This operator-() const {
    This ret;
    ret.iv_= -iv_;
    if (rep_) {
      ret.rep_= new Rep(-rep_->isolated_root());
    }
  }
  double to_double() const {
    if (!rep_) return iv_.first;
    else return to_double(rep_->isolated_root());
  }
  std::pair<double,double> to_interval() const {
    if (!rep_) return iv_;
    else return to_interval(rep_->isolated_root());
  }
  bool is_rational() const {
    return rep_->isolated_root().is_rational();
  }
  static This infinity() {
    This ret;
    ret.iv_= std::numeric_limits<double>::infinity();
  }
protected:
  Interval_nt iv_;
  typename Rep::Shared_root root_;
  typename Rep::Pointer rep_;  
};

template <class Traits_t>
class Lazy_upper_bound_root_stack{
  typedef typename Traits_t::Root TRoot;

  struct Rep: public Ref_counted<Rep> {
    
    TRoot root_;
  };

 
public:
  typedef Traits_t Traits;
  typedef Lazy_upper_bound_root_stack_root<Traits, Rep> Root;
  typedef typename Traits::Function Function;

  Upper_bound_root_stack(): cur_(Root::infinity()) {
  };
  Upper_bound_root_stack(const Polynomial &f, 
			 const Root &lb,
			 const Root &ub,
			 const Traits &tr): rep_(new Rep(f, lb, ub, tr)){}
  
  void pop() {
    rep_->pop();
  }
  const Root &top() const {
    return rep_->top();
  }

  bool empty() const {
    return rep_->empty();
  }

protected:
  Rep::Pointer rep_;
};

CGAL_POLYNOMIAL_END_NAMESPACE


#endif
