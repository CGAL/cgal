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

#include <CGAL/Polynomial/internal/numeric_solvers_support.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>

#ifdef CGAL_HAVE_TNT
#include <TNT/tnt_array2d.h>
#include <TNT/tnt_array1d.h>
#include <TNT/jama_eig.h>
#endif

#include <algorithm>
#include <functional>
#include <iterator>

//#include <iomanip>

namespace CGAL { namespace POLYNOMIAL { namespace internal {
#if CGAL_HAVE_TNT
//static const double max_error_value =0.00005;

template <bool CLEAN, class NT>
static void jama_compute_roots(const NT *begin, const NT *end,  NT lb,
NT ub, std::vector<NT> &roots)
{
    int degree= end-begin-1;
    TNT::Array2D<NT> arr(degree, degree, 0.0);
    for (int i=0; i< degree; ++i) {
        arr[0][i]=-begin[degree-i-1]/begin[degree];
    }
    for (int i=0; i+1< degree; ++i) {
        arr[i+1][i]=1;
    }

    JAMA::Eigenvalue<NT> ev(arr);
    TNT::Array1D<NT> real, imag;
    ev.getImagEigenvalues(imag);
    ev.getRealEigenvalues(real);
    CGAL_Polynomial_assertion(imag.dim1()== real.dim1());

    /*NT tol;
    if (CLEAN) tol=.00005;
    else tol=0;*/

    for (int i=0; i< real.dim1(); ++i) {
        if (root_is_good(real[i], imag[i], lb-tol, ub)) {
            roots.push_back(real[i]/*polish_root(begin, end, real[i])*/);
        } else {
        }
    }
    std::sort(roots.begin(), roots.end(), std::greater<NT>());
    if (CLEAN) filter_roots(begin, end, lb, roots);
}

#endif



void jama_polynomial_compute_roots(const double *begin, const double *end,
double lb, double ub,
std::vector<double> &roots)
{
  std::ptrdiff_t degree= end-begin-1;
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
        default:
#ifdef CGAL_HAVE_TNT
	  jama_compute_roots<false>(begin, end, lb, ub, roots);
#else
	  CGAL_error();
#endif
	  //jama_compute_roots<false>(begin, end, lb, ub, roots);
    }
}


void jama_polynomial_compute_cleaned_roots(const double *begin, const double *end,
double lb, double ub,
std::vector<double> &roots)
{
  std::ptrdiff_t degree= end-begin-1;
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
        default:
#ifdef CGAL_HAVE_TNT
	  jama_compute_roots<true>(begin, end, lb, ub, roots);
#else
	  CGAL_error();
#endif
    }
}


} } } //namespace CGAL::POLYNOMIAL::internal
