// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_NUMERIC_ROOT_STACK_H
#define CGAL_POLYNOMIAL_NUMERIC_ROOT_STACK_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/internal/Double_with_infinity.h>
#include <iterator>
#ifdef CGAL_USE_GSL
#include <CGAL/Polynomial/internal/GSL_numeric_solver.h>
#endif

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <bool HINF=std::numeric_limits<double>::has_infinity>
  struct Numeric_pick_root {
    typedef double Root;
    //typedef internal::Double_with_infinity Root;
  };

  template <>
  struct Numeric_pick_root<false> {
    typedef internal::Double_with_infinity Root;
  };
} } } //namespace CGAL::POLYNOMIAL::internal

namespace CGAL { namespace POLYNOMIAL {

#ifdef CGAL_USE_GSL
#define CGAL_DEFAULT_NUMERIC_SOLVER CGAL::POLYNOMIAL::internal::GSL_numeric_solver
#define CGAL_DEFAULT_CLEANED_NUMERIC_SOLVER CGAL::POLYNOMIAL::internal::GSL_cleaned_numeric_solver
#else
#define CGAL_DEFAULT_NUMERIC_SOLVER CGAL::POLYNOMIAL::internal::Turkowski_numeric_solver
#define CGAL_DEFAULT_CLEANED_NUMERIC_SOLVER CGAL::POLYNOMIAL::internal::Turkowski_cleaned_numeric_solver
#endif

template <class Solver_traits, class Numeric_solver=CGAL_DEFAULT_NUMERIC_SOLVER >
class Numeric_root_stack
{
public:
    typedef internal::Numeric_pick_root<>::Root Root;
protected:

  /* All this mess is to handle when roots are not doubles and when the coefficients are not doubles
   */

  template <class Rt>
  void initialize_2(const double *b, const double *e, double lb, double ub, std::vector<Rt> &roots) {
    std::vector<double> lroots;
    ns_(b, e, lb, ub, lroots);
    roots.insert(roots.end(), lroots.begin(), lroots.end());
  }

  void initialize_2(const double *b, const double *e, double lb, double ub, std::vector<double> &roots) {
    ns_(b, e, lb, ub, roots);
  }

  template <class Fn>
  void initialize(const Fn &f, Root lb, Root ub) {
    std::vector<double> c(f.degree()+1);
    for (unsigned int i=0; i<= f.degree(); ++i) {
      c[i]= to_double(f[i]);
    }
    initialize_2(&*c.begin(), &*c.begin()+ f.degree()+1, 
		 static_cast<double>(lb), static_cast<double>(ub),
		 roots_);
    /*if (CLEAN) {
      polynomial_compute_cleaned_roots(&*c.begin(), &*c.begin()+ f.degree()+1, lb, ub, roots_);
      } else {
      polynomial_compute_roots(&*c.begin(), &*c.begin()+ f.degree()+1, lb, ub, roots_);
      }*/
  }
  void initialize(const Polynomial<double> &f, Root lb, Root ub) {
    const double *p0= &*f.begin();
    initialize_2(p0, p0+ f.degree()+1, static_cast<double>(lb), static_cast<double>(ub), roots_);
    /*if (CLEAN) {
      polynomial_compute_cleaned_roots(&*f.begin(), &*f.begin()+ f.degree()+1, lb, ub, roots_);
      } else {
      polynomial_compute_roots(&*f.begin(), &*f.begin()+ f.degree()+1, lb, ub, roots_);
      }*/
  }


public:

  typedef typename Solver_traits::Function Function;
  typedef Solver_traits Traits;
  Numeric_root_stack(const typename Solver_traits::Function &f, Root lb, Root ub, const Solver_traits&) {
    //std::cout << "Solving " << f << " from " << lb << " to " << ub;
    initialize(f, lb, ub);	 
    /*if (!roots_.empty()) std::cout << " got " << roots_.back() << std::endl;
      else std::cout << std::endl;*/
#if 0
    for (unsigned int i=1; i < roots_.size();  ++i) {
      if (roots_[i]>  roots_[i-1]){
	std::cerr << "ERROR: roots out of order ";
	std::copy(roots_.begin(), roots_.end(), std::ostream_iterator<double>(std::cerr, " "));
	std::cerr << " for " << f << " from " << lb << " to " << ub << std::endl;
	roots_.clear();
	initialize(f,lb, ub);
      }
    }
#endif
    //CGAL_Polynomial_postcondition(roots_[i] <= roots_[i-1]);
  }

  Numeric_root_stack(){};

  void pop() {
    CGAL_Polynomial_precondition(!roots_.empty());
    roots_.pop_back();
  }

  const Root& top() const
  {
    CGAL_Polynomial_precondition(!roots_.empty());
    return roots_.back();
  }

  bool empty() const
  {
    return roots_.empty();
  }

protected:
  Numeric_solver ns_;
  std::vector<Root> roots_;
};

} } //namespace CGAL::POLYNOMIAL
#endif
