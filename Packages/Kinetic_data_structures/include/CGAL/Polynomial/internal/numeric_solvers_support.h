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

#ifndef CGAL_POLYNOMIAL_INTERNAL_NUMERIC_SOLVERS_SUPPORT_H
#define CGAL_POLYNOMIAL_INTERNAL_NUMERIC_SOLVERS_SUPPORT_H
#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

void compute_quadratic_roots(const double *begin, const double *end, 
			     double lb, double ub, std::vector<double> &roots);

void compute_quadratic_cleaned_roots(const double *begin, const double *end, 
				     double lb, double ub, std::vector<double> &roots);

void compute_linear_roots(const double *begin, const double *end, 
			  double lb, double ub, std::vector<double> &roots);

void compute_linear_cleaned_roots(const double *begin, const double *end, 
				  double lb, double ub, std::vector<double> &roots);

void check_first_root(const double *begin, const double *end, 
		      double lb, std::vector<double> &roots);




template <class NT>
inline bool root_is_good(NT r, NT c, NT lb, NT ub, NT tol=.0000005){
  if (std::abs(c) < tol && r > lb && r < ub) {
    return true;
  } else return false;
}

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
