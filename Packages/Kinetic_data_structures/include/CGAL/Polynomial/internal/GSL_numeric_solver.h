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

#ifndef CGAL_POLYNOMIAL_INTERNAL_GSL_NUMERIC_SOLVER_H
#define CGAL_POLYNOMIAL_INTERNAL_GSL_NUMERIC_SOLVER_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/numeric_solvers_support.h>
#include <vector>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <bool CLEAN>
inline double gsl_max_error() {
  if (CLEAN) return .0000005;
  else return 0;
}

template <bool CLEAN>
inline void gsl_compute_roots(const double *begin, const double *end, 
				     double lb, double ub, 
				     std::vector<double> &roots_){

  double max_error=gsl_max_error<CLEAN>();

  gsl_poly_complex_workspace *workspace;
  int degree= end-begin-1;
  workspace= gsl_poly_complex_workspace_alloc(degree+1);
    
  //! I can't use auto_ptr because it is an array
  double *roots= new double[2*degree];
    
  /*for ( int i=0; i< 2*fn.degree(); ++i){
    roots[i]=0;
    }*/
    
  int ret = gsl_poly_complex_solve(begin, degree+1, workspace,
				   roots);
  roots_.reserve(ret);

  if (ret!=GSL_SUCCESS){
    std::cerr << "GSL solver did not converge on ";
    std::copy(begin, end, std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << std::endl;
  }
    
  for ( int i=degree-1; i>=0; --i){
    double r= roots[2*i];
    double c= roots[2*i+1];
    if (root_is_good(r, c, lb-max_error, ub)) {
      roots_.push_back(r);
    }
  }

  std::sort(roots_.begin(), roots_.end(), std::greater<double>() );
  delete roots;

  gsl_poly_complex_workspace_free(workspace);

  if (CLEAN) check_first_root(begin,end, lb, roots_);
  return;
}





template <bool CLEAN>
inline void gsl_compute_cubic_roots(const double *begin, const double *end, 
				    double lb, double ub, 
				    std::vector<double> &roots_){
  CGAL_Polynomial_precondition(begin[3] != 0);
  double max_error=gsl_max_error<CLEAN>();

  double r[3];
  int num_roots= gsl_poly_solve_cubic(begin[2]/begin[3], 
				      begin[1]/begin[3],
				      begin[0]/begin[3], &r[0],&r[1],&r[2]);
  roots_.reserve(num_roots);
  // I want reverse sorted roots
  for (int i=num_roots-1; i>=0; --i){
    if (r[i]>  lb- max_error && r[i] < ub) roots_.push_back(r[i]);
  }
  if (CLEAN) check_first_root(begin, end, lb, roots_);
}


inline void gsl_polynomial_compute_roots(const double *begin, const double *end, 
				  double lb, double ub, 
				  std::vector<double> &roots){
   int degree= end-begin-1;
  switch( degree) {
  case -1:
  case 0:
    break;
  case 1:
    compute_linear_roots(begin,end, lb, ub, roots);
    break;
  case 2:
    compute_quadratic_roots(begin, end, lb, ub, roots);
    break;
  case 3:
    gsl_compute_cubic_roots<false>(begin, end, lb, ub, roots);
    break;
  default:
    gsl_compute_roots<false>(begin, end, lb, ub, roots);
  }
}

inline void gsl_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
				  double lb, double ub, 
				  std::vector<double> &roots){
  int degree= end-begin-1;
  switch( degree) {
  case -1:
  case 0:
    break;
  case 1:
    compute_linear_cleaned_roots(begin,end, lb, ub, roots);
    break;
  case 2:
    compute_quadratic_cleaned_roots(begin, end, lb, ub, roots);
    break;
  case 3:
    gsl_compute_cubic_roots<true>(begin, end, lb, ub, roots);
    break;
  default:
    gsl_compute_roots<true>(begin, end, lb, ub, roots);
  }
}

inline double gsl_evaluate_polynomial(const double *b, const double *e, double t) {
  if (b==e) return 0;
  /*double *d= new double[coefs.size()];
  for (unsigned int i=0; i< coefs.size(); ++i){
    d[i]= coefs[i];
    }*/
  double v= gsl_poly_eval(b, e-b, t);
  return v;
}
/*

inline void gsl_polynomial_compute_roots(const double *begin, const double *end, 
				  double lb, double ub, 
				  std::vector<double> &roots){
 CGAL_assertion(0&& no_gsl_support_built);
}

inline void gsl_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
				  double lb, double ub, 
				  std::vector<double> &roots){
  CGAL_assertion(0&& no_gsl_support_built);
}

inline double gsl_evaluate_polynomial(const double *b, const double *e, double t) {
  CGAL_assertion(0&& no_gsl_support_built);
  return std::numeric_limits<double>::nan();
}
*/

/*// GSL
void gsl_polynomial_compute_roots(const double *begin, const double *end, 
				  double lb, double ub, std::vector<double> &roots);

void gsl_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
double lb, double ub, std::vector<double> &roots);*/



struct GSL_numeric_solver {
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
  };


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
