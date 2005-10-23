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

#ifndef CGAL_POLYNOMIAL_NUMERIC_ROOT_STACK_H
#define CGAL_POLYNOMIAL_NUMERIC_ROOT_STACK_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

template <class Solver_traits, class Numeric_solver=internal::Turkowski_numeric_solver>
class Numeric_root_stack {

  template <class Fn>
  void initialize(const Fn &f, double lb, double ub) {
    std::vector<double> c(f.degree()+1);
    for (unsigned int i=0; i<= f.degree(); ++i){
      c[i]= to_double(f[i]);
    }
    ns_(&*c.begin(), &*c.begin()+ f.degree()+1, lb, ub, roots_);
    /*if (CLEAN) {
      polynomial_compute_cleaned_roots(&*c.begin(), &*c.begin()+ f.degree()+1, lb, ub, roots_);
    } else {
      polynomial_compute_roots(&*c.begin(), &*c.begin()+ f.degree()+1, lb, ub, roots_);
      }*/
  }
  void initialize(const Polynomial<double> &f, double lb, double ub) {
    const double *p0= &*f.begin();
    ns_(p0, p0+ f.degree()+1, lb, ub, roots_);
    /*if (CLEAN) {
      polynomial_compute_cleaned_roots(&*f.begin(), &*f.begin()+ f.degree()+1, lb, ub, roots_);
    } else {
      polynomial_compute_roots(&*f.begin(), &*f.begin()+ f.degree()+1, lb, ub, roots_);
      }*/
  }
public:
  typedef double Root;
  typedef typename Solver_traits::Function Function;
  typedef Solver_traits Traits;
  Numeric_root_stack(const typename Solver_traits::Function &f, Root lb, Root ub, const Solver_traits&){
    initialize(f, lb, ub);
    for (unsigned int i=1; i < roots_.size();  ++i){
      CGAL_Polynomial_postcondition(roots_[i] <= roots_[i-1]);
    }
  }
  
  Numeric_root_stack(){};

  void pop() {
    CGAL_Polynomial_precondition(!roots_.empty());
    roots_.pop_back();
  }

  const Root& top() const {
    CGAL_Polynomial_precondition(!roots_.empty());
    return roots_.back();
  }

  bool empty() const {
    return roots_.empty();
  }

protected:
  Numeric_solver ns_;
  std::vector<double> roots_;
};

CGAL_POLYNOMIAL_END_NAMESPACE
#endif




