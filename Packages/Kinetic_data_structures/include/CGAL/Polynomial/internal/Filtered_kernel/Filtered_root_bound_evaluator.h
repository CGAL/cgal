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

#ifndef CGAL_POLYNOMIAL_FILTERED_ROOT_BOUND_EVALUATOR_H
#define CGAL_POLYNOMIAL_FILTERED_ROOT_BOUND_EVALUATOR_H

#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE


template<class Kernel, class M_t = CGAL::Field_tag>
class Filtered_root_bound_evaluator
{
public:
  Filtered_root_bound_evaluator(bool pow,
				const Kernel k): rb_(k.interval_kernel_object().root_bound_object(pow))  {}
  typedef double result_type;
  typedef typename Kernel::Function argument_type;
  
  result_type operator()(const argument_type& p) const
  {
    Interval_arithmetic_guard iag;
    return rb_(p.interval_function()).sup();
  }
protected:
  typename Kernel::Interval_kernel::Root_bound rb_;
};



CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif // CGAL_POLYNOMIAL_ROOT_BOUND_EVALUATOR_H
