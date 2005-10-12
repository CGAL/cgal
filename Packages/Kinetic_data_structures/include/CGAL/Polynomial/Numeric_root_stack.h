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




