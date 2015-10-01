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

#ifndef CGAL_POLYNOMIAL_INTERNAL_NUMERIC_SOLVERS_SUPPORT_H
#define CGAL_POLYNOMIAL_INTERNAL_NUMERIC_SOLVERS_SUPPORT_H
#include <CGAL/Polynomial/basic.h>
#include <vector>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class NT>
inline bool root_is_good(NT r, NT c, NT lb, NT ub, NT tol=.0000005)
{
  if (std::abs(c) < tol && r > lb && r < ub) {
        return true;
    } else return false;
}


void compute_quadratic_roots(const double *begin, const double *end,
			     double lb, double ub,
			     std::vector<double> &roots);

void compute_quadratic_cleaned_roots(const double *begin, const double *end,
				     double lb, double ub, 
				     std::vector<double> &roots);

void compute_linear_roots(const double *begin, const double *end,
			  double lb, double ub,
			  std::vector<double> &roots);

void compute_linear_cleaned_roots(const double *begin, const double *end,
				  double lb, double ub, 
				  std::vector<double> &roots);



void filter_solver_roots(const double *begin, const double *end,
			 double lb, double ub, double last,
			 std::vector<double> &roots);

} } } //namespace CGAL::POLYNOMIAL::internal

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Polynomial/internal/numeric_solvers_support_impl.h>
#endif // CGAL_HEADER_ONLY

#endif
