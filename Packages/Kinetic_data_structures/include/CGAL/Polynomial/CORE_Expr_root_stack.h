#ifndef POLYNOMIAL_CORE_SOLVER_H
#define POLYNOMIAL_CORE_SOLVER_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Polynomial/internal/Explicit_root.h>
#include <CGAL/NT_converter.h>
#include <CORE/BigInt.h>

#include <iostream>

CGAL_BEGIN_NAMESPACE 
double to_double(const CORE::BigInt &bi){
  return bi.doubleValue();
}
CGAL_END_NAMESPACE

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

template <class Solver_traits>
class CORE_Expr_root_stack {
protected:
  typedef typename Solver_traits::NT Coef;
  typedef typename CORE::Polynomial<Coef> CORE_polynomial;
  typedef typename CORE::Sturm<Coef> CORE_Sturm;
  typedef CORE_Expr_root_stack<Solver_traits> This;
public:
  typedef Solver_traits Traits;
  typedef Explicit_root<CORE::Expr> Root;

  //! NOTE: The function must be square free!!!!!!!!!!!
  CORE_Expr_root_stack(const typename Solver_traits::Function &f,  
		       const Root &lb, 
		       const Root &ub,
		       const Traits &): ub_(ub) {
    initialize(f, lb);
  }

  CORE_Expr_root_stack(): counter_(1), num_roots_(0){}

  const Root& top() const {
    return cur_;
  }
  void pop() {
    Polynomial_precondition(counter_<=num_roots_);
    if (counter_ == num_roots_) {
      no_roots();
    } else {
      ++counter_;
      make_cur_root(make_expr());
      enforce_upper_bound();
    }
  }

  bool empty() const {
    return counter_>num_roots_;
  }
  
protected:
  CORE_polynomial poly_;
  Root  ub_;
  Root cur_;
  int counter_;
  int num_roots_;
  
  void initialize(const typename Solver_traits::Function &fn, const Root& lb){
    if (fn.degree()<0){
      no_roots();
      return;
    } else {
      poly_= CORE_polynomial(fn.degree());
      for (int i=0; i<= fn.degree(); ++i){
	poly_.setCoeff(i, fn[i]);
      }
      initialize_counters(lb);
	
      CORE::Expr testr;
      do {
	++counter_;
	if (counter_> num_roots_){
	  no_roots();
	  return;
	}
	testr= make_expr();
	
      } while (lb != -infinity<Root>() && testr <= lb.representation());
      make_cur_root(testr);
      
    }
    //std::cout << "There are " << _num_roots << " roots.\n";
    //std::cout << "Counter is set to " << _counter << "\n";
    enforce_upper_bound();
  }

  void enforce_upper_bound(){
    if (cur_ < ub_) return;
    else no_roots();
  }

  CORE::Expr make_expr() {
    return CORE::Expr(CORE::rootOf(poly_, counter_));
  }
  void make_cur_root(const CORE::Expr &r){
    cur_= Root(r);
  }
  
  void no_roots(){
    ub_= CORE::Expr(0);
    cur_= infinity<Root>();
    num_roots_=0;
    counter_=1;
  }

  void initialize_counters(const Root &lb){
    std::cout << "Computing strum of " << poly_ << "..." << std::flush;
    CORE_Sturm sturm(poly_);
    std::cout << "done." << std::endl;
    num_roots_=0;
    CGAL_assertion(-ub_ != infinity<Root>());
    if (lb== -infinity<Root>() && ub_== infinity<Root>()){
      num_roots_= sturm.numberOfRoots();
      counter_=0;
    } else if (ub_ == infinity<Root>()){
      num_roots_= sturm.numberOfRootsAbove(bf_lower_bound(lb.representation()));
      counter_ = sturm.numberOfRootsBelow(bf_lower_bound(lb.representation()));
    } else if (lb == infinity<Root>()){
      num_roots_= sturm.numberOfRootsBelow(bf_upper_bound(ub_.representation()));
      counter_ = 0;
    } else {
      counter_= sturm.numberOfRootsBelow(bf_lower_bound(lb.representation()));
      num_roots_= sturm.numberOfRoots(bf_lower_bound(lb.representation()),
				      bf_upper_bound(ub_.representation()));
    }
  };

  //! There are probably better ways of doing this
  Coef bf_lower_bound(const CORE::Expr &rt) const {
    machine_double lb, ub;
    rt.doubleInterval(lb, ub);
    return Coef(lb);
  }

  //! There are probably better ways of doing this
  Coef bf_upper_bound(const CORE::Expr &rt) const {
    machine_double lb, ub;
    rt.doubleInterval(lb, ub);
    return Coef(ub);
  }
};

CGAL_POLYNOMIAL_END_NAMESPACE;

#endif
