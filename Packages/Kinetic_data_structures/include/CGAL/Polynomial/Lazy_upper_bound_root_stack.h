// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_LAZY_UPPER_BOUND_ROOT_STACK_H
#define CGAL_POLYNOMIAL_LAZY_UPPER_BOUND_ROOT_STACK_H
#include <CGAL/Polynomial/Upper_bound_root_stack.h>
#include <CGAL/KDS/Ref_counted.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

/*
  Root has pointer to Root_rep, solver_rep (ref counted), a double approximation;

  Root_rep is 
  - an optional Root
  - a pointer to the solver (non-ref_counted)

  Solver_rep is
  - an normal solver
  - pointer to current root (ref counted)
  
  Top does:
  - creates a new root rep with a pointer to the current root

  Pop does:
  - isolates current root if necessary
  - puts it in to Root_rep
  - creates new Root_rep (if more roots)

  Compare:
  - compares double approximations
  - if root is not isolated, isolates root,
  - updates double approx
  - compares double approximation
  - compares isolated root


  Optimizations: 

  - keep popping intervals off of the front of the solver to try to
  perform the comparison, stop when isolated. Try this only initially with a nt or an isolated root.
*/


template <class Rep, class Solver_rep>
class Lazy_upper_bound_root_stack_root {
  typedef Lazy_upper_bound_root_stack_root<Rep, Solver_rep> This;

  Comparison_result compare_double(const This &o) const {
    if (iv_.inf() > o.iv_.sup()) return LARGER;
    else if (iv_.sup() < o.iv_.inf()) return SMALLER;
    else if (iv_.is_point() && o.iv_.is_point()) return EQUAL;
    else return UNKNOWN;
  }

  Comparison_result compare(const This &o) const {
    int cd= compare_double(o);
    if (cd != UNKNOWN) {
      return cd;
    } else {
      update_interval();
      o.update_interval();
      
      int cd= compare_double(o);
      if (cd != UNKNOWN){
	return cd;
      } else {
	if (rep_->has_root() && !o.rep_->has_root()) {
	  if (o.srep_->refine_top(iv_.sup())) {
	    return SMALLER;
	  }
	} else if (o.rep_->has_root() && !rep_->has_root()) {
	  if (srep_->refine_top(o.iv_.sup())) {
	    return LARGER;
	  }
	}
	ensure_exact();
	o.ensure_exact();
	return rep_->root().compare(o.rep_->root());
      }
    }
    CGAL_Polynomial_postcondition(0);
    return  UNKNOWN;
  }

  
public:
  template <class NT>
  Lazy_upper_bound_root_stack_root(const NT &nt): iv_(CGAL_POLYNOMIAL_TO_INTERVAL(nt)), rep_(new Rep(nt)) {}
  Lazy_upper_bound_root_stack_root(double d): iv_(d), rep_(new Rep(d)) {}

  Lazy_upper_bound_root_stack_root(typename Solver_rep::Pointer sp): iv_(sp->double_interval()), srep_(sp), 
							    rep_(sp->current_root_rep()){}
  Lazy_upper_bound_root_stack_root(){}
  
  bool operator==(const This &o) const {
    return compare(o)==EQUAL;
  }
  bool operator<(const This &o) const {
    return compare(o)==SMALLER;
  }
  bool operator>(const This &o) const {
    return compare(o)==LARGER;
  }
  bool operator>=(const This &o) const {
    return compare(o)!=SMALLER;
  }
  bool operator<=(const This &o) const {
    return compare(o)!=LARGER;
  }
  const typename Rep::Root &troot() const {
    ensure_exact();
    return rep_->root();
  }
  bool is_even_multiplicity() const {
    ensure_exact();
    return rep_->root().is_even_multiplicity();
  }
  This operator-() const {
    ensure_exact();
    This ret;
    ret.iv_= -iv_;
    ret.rep_= new Rep(-rep_->root());
    return ret;
  }
  double to_double() const {
    ensure_exact();
    return CGAL_POLYNOMIAL_TO_DOUBLE(rep_->root());
  }
  std::pair<double,double> double_interval(double tol= std::numeric_limits<double>::infinity()) const {
    if (tol != std::numeric_limits<double>::infinity()) {
      ensure_exact();
      iv_= CGAL_POLYNOMIAL_TO_INTERVAL(rep_->root());
    } else {
      return iv_;
    }
  }
  bool is_rational() const {
    return rep_->root(srep_).is_rational();
  }
  static This infinity() {
    This ret;
    ret.iv_= std::numeric_limits<double>::infinity();
  }
  
  void write(std::ostream &out) const {
    if (!rep_) out << "null";
    if (rep_->has_root()) out << rep_->root();
    else out << "Unisolated (" << iv_ << ")";
  }

protected:
  void ensure_exact() const {
    if (!rep_->has_root()){
      srep_->check_current(rep_);
      srep_->isolate();
    }
    iv_= CGAL_POLYNOMIAL_TO_INTERVAL(rep_->root());
  }

  void update_interval() const {
    if (rep_->has_root()){
      iv_= rep_->root().double_interval(std::numeric_limits<double>::infinity());
    } else {
      srep_->check_current(rep_);
      iv_= srep_->double_interval();
    }
  }

  mutable Interval_nt iv_;
  mutable typename Solver_rep::Pointer srep_;  
  mutable typename Rep::Pointer rep_;
};


template <class Rep, class Solver_rep>
std::ostream& operator<<(std::ostream &out, const Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt){
  rt.write(out);
  return out;
}

template <class Rep, class Solver_rep>
double to_double(const Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt){
  return rt.to_double();
}

template <class Rep, class Solver_rep>
std::pair<double,double> to_interval(const Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt){
  return rt.double_interval();
}


template <class Traits_t>
class Lazy_upper_bound_root_stack{
  typedef typename Traits_t::Root TRoot;
protected:
  typedef Upper_bound_root_stack<Traits_t> Solver;
  struct Solver_rep;
  struct Root_rep;
public:
  typedef Traits_t Traits;
  typedef Lazy_upper_bound_root_stack_root<Root_rep, Solver_rep> Root;
  typedef typename Traits::Function Function;

  Lazy_upper_bound_root_stack(){
  };
  Lazy_upper_bound_root_stack(const Function &f, 
			      const Root &lb,
			      const Root &ub,
			      const Traits &tr): rep_(new Solver_rep(f, lb, ub, tr)){}
  
  void pop() {
    rep_->pop();
  }
  Root top() const {
    return rep_;
  }

  bool empty() const {
    return rep_->empty();
  }

protected:
  mutable typename Solver_rep::Pointer rep_;
};

template <class Traits_t>
struct Lazy_upper_bound_root_stack<Traits_t>::Solver_rep: public CGAL::KDS::Ref_counted<Solver_rep> {
  Solver_rep(const Function &f, 
	     const Root &lb,
	     const Root &ub,
	     const Traits &tr): solver_(f, lb.troot(), ub.troot(), tr, false){
    cur_= new Root_rep();
  }

  void check_current(typename Root_rep::Pointer cur) const {
    CGAL_assertion(cur== cur_);
  }

  typename Root_rep::Pointer current_root_rep() const {   
    CGAL_precondition(cur_);
    return cur_;
  }

  typename Root_rep::Pointer current_root_rep() {
    CGAL_precondition(cur_);
    return cur_;
  }
  
  void pop() {
    if (!is_isolated()) isolate();
    
    if (!solver_.empty()) {
      cur_= new Root_rep();
    } else {
      cur_=NULL;
    }
  }

  std::pair<double, double> double_interval() const {
    return solver_.double_interval();
  }

  bool empty() const {
    return !cur_;
  }

  void isolate() {
    solver_.isolate_top();
    set_root();
  }

  bool is_isolated() const {
    return cur_->has_root();
  }

  bool refine_bottom(double ub) {
    bool ret= solver_.refine_bottom(ub);
    if (solver_.top_is_isolated()){
      set_root();
    }
    return ret;
  }

  bool refine_top(double lb) {
    bool ret= solver_.refine_top(lb);
    if (solver_.top_is_isolated()){
      set_root();
    }
    return ret;
  }

protected:
  void set_root() {
    CGAL_precondition(!cur_->has_root());
    cur_->set_root(solver_.top());
    solver_.pop_no_isolate();
  }

  Solver solver_;
  typename Root_rep::Pointer cur_;
};

template <class Traits_t>
struct Lazy_upper_bound_root_stack<Traits_t>::Root_rep: public CGAL::KDS::Ref_counted<Root_rep> {
  typedef TRoot Root;
  Root_rep(): root_(-std::numeric_limits<TRoot>::infinity()) {
  }
  Root_rep(const TRoot &tr): root_(tr){}
  template <class NT>
  Root_rep(const NT &nt): root_(nt){}

  Root_rep(double d): root_(d){}

  void set_root(const TRoot &tr) {
    CGAL_precondition(root_==-std::numeric_limits<TRoot>::infinity());
    root_=tr;
  }

  const TRoot& root() const {
    CGAL_precondition(root_!= -std::numeric_limits<TRoot>::infinity());
    return root_;
  }
  
  bool has_root() const {
    return root_!= -std::numeric_limits<TRoot>::infinity();
  }

  /*std::pair<double,double> current_interval(typename Solver_rep::Pointer srep) const {
    if (is_isolated()){
      return root_.double_interval();
    } else {
      srep->check_current(this);
      return srep->double_interval();
    }
    }*/
  
  mutable TRoot root_;
};

CGAL_POLYNOMIAL_END_NAMESPACE

CGAL_BEGIN_NAMESPACE
template <class Rep, class Solver_rep>
double to_double(const CGAL_POLYNOMIAL_NS::Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt){
  return rt.to_double();
}

template <class Rep, class Solver_rep>
std::pair<double,double> to_interval(const CGAL_POLYNOMIAL_NS::Lazy_upper_bound_root_stack_root<Rep, Solver_rep> &rt){
  return rt.double_interval();
}
CGAL_END_NAMESPACE

#endif
