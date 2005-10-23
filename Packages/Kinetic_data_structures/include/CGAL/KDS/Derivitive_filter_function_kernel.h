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

#ifndef CGAL_KDS_NUMERIC_SOLVER_H
#define CGAL_KDS_NUMERIC_SOLVER_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Numeric_root_stack.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>


CGAL_KDS_BEGIN_NAMESPACE



template <class Traits>
struct Derivitive_filter_function_kernel: public Traits {
  
  typedef CGAL::POLYNOMIAL::Numeric_root_stack<Traits,
					       CGAL::POLYNOMIAL::internal::Turkowski_cleaned_numeric_solver> Root_stack;
  /*class Root_stack: public CGAL_POLYNOMIAL_NS::internal::Numeric_root_stack_core<Traits, true> {
    typedef CGAL_POLYNOMIAL_NS::internal::Numeric_root_stack_core<Traits, true> Parent;
  public:
    typedef typename Parent::Root Root;
    typedef typename Traits::Function Function;
    Root_stack(const typename Traits::Function &f, 
	       Root lb, Root ub, 
	       const Traits&k): Parent(f, lb, ub, k){
      CGAL_KDS_LOG(LOG_LOTS, "Solved " << f << " from " << lb << " to " << ub << " to get ");
      for (unsigned int i=0; i< Parent::roots_.size(); ++i){
	CGAL_KDS_LOG(LOG_LOTS, Parent::roots_[i] << " ");
      }
      CGAL_KDS_LOG(LOG_LOTS, std::endl);
    }

    Root_stack(){};
    };*/

  typedef typename Root_stack::Root Root;

  Derivitive_filter_function_kernel(Traits tr): Traits(tr){}
  Derivitive_filter_function_kernel(){}

  Root_stack root_stack_object(const typename Traits::Function &f,
			       const Root &lb,
			       const Root &ub) {
    return Root_stack(f, lb, ub, *this);
  }

};

CGAL_KDS_END_NAMESPACE




#endif // inclusion guard
