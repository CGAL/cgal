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

#ifndef CGAL_POLYNOMIAL_INTERNAL_NUMERIC_SOLVERS_H
#define CGAL_POLYNOMIAL_INTERNAL_NUMERIC_SOLVERS_H
#include <CGAL/Polynomial/basic.h>
#include <vector>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

// JAMA
void jama_polynomial_compute_roots(const double *begin, const double *end, 
				   double lb, double ub, std::vector<double> &roots);

void jama_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
					   double lb, double ub, std::vector<double> &roots);

// GSL
void gsl_polynomial_compute_roots(const double *begin, const double *end, 
				  double lb, double ub, std::vector<double> &roots);

void gsl_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
					  double lb, double ub, std::vector<double> &roots);

// Turkowski
void Turkowski_polynomial_compute_roots(const double *begin, const double *end, 
				  double lb, double ub, std::vector<double> &roots);

void Turkowski_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
					  double lb, double ub, std::vector<double> &roots);

/*struct GSL_numeric_solver {
  void operator()(const double *begin, const double *end,  
		  double lb, double ub,
		  std::vector<double> &roots) const {
    gsl_polynomial_compute_roots(begin, end, lb, ub, roots);
  }
};

struct GSL_cleaned_numeric_solver {
  void operator()(const double *begin, const double *end,  
		  double lb, double ub,
		  std::vector<double> &roots) const {
    gsl_polynomial_compute_cleaned_roots(begin, end, lb, ub, roots);
  }
  };*/

struct JAMA_numeric_solver {
  void operator()(const double *begin, const double *end,  
		  double lb, double ub,
		  std::vector<double> &roots) const {
    jama_polynomial_compute_roots(begin, end, lb, ub, roots);
  }
};

struct JAMA_cleaned_numeric_solver {
  void operator()(const double *begin, const double *end,  
		  double lb, double ub,
		  std::vector<double> &roots) const {
    jama_polynomial_compute_cleaned_roots(begin, end, lb, ub, roots);
  }
};

struct Turkowski_numeric_solver {
  void operator()(const double *begin, const double *end,  
		  double lb, double ub,
		  std::vector<double> &roots) const {
    Turkowski_polynomial_compute_roots(begin, end, lb, ub, roots);
  }
};

struct Turkowski_cleaned_numeric_solver {
  void operator()(const double *begin, const double *end,  
		  double lb, double ub,
		  std::vector<double> &roots) const {
    Turkowski_polynomial_compute_cleaned_roots(begin, end, lb, ub, roots);
  }
};



CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
