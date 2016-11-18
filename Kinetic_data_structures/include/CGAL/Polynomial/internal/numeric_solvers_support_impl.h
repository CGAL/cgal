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


#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/Polynomial/internal/numeric_solvers_support.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>

#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/internal/Rational/Derivative.h>
#include <CGAL/Polynomial/Interval_polynomial.h>

/*#ifdef _MSC_VER
#pragma warning(disable:1572)
#endif*/

namespace CGAL { namespace POLYNOMIAL { namespace internal {

static const double max_error_value=.00005;

namespace {
template <bool CLEAN, class NT>
inline void compute_quadratic_roots_t(const NT *begin, const NT * /*end*/,  NT lb, NT ub,
std::vector<NT> &roots)
{
  NT max_error=0;
  if (CLEAN) max_error=max_error_value;
  CGAL_Polynomial_assertion(begin[2] != 0);
  
  NT desc= begin[1]*begin[1]-4*begin[0]*begin[2];
  if (desc <= 0) return;
  
  NT ur= (-begin[1]+sqrt(desc))/(2*begin[2]);
  NT lr= (-begin[1]-sqrt(desc))/(2*begin[2]);
  if (begin[2]< 0) std::swap(lr, ur);
  if (lr > lb-max_error && lr < ub) {
    roots.push_back(ur);
    if (lr > lb-max_error && lr < ub && (!CLEAN || /*lr > lb+max_error ||*/ begin[2] >0)){
      roots.push_back(lr);
    }
  } else {
    // only upper
    if (ur > lb-max_error && ur < ub && (!CLEAN || /*ur > lb+max_error ||*/ begin[2] <0)){
      roots.push_back(ur);
    }


  }

  // drop even roots
  /*if (ur >lb-max_error && ur < ub){
    if (!CLEAN || sign(begin[2]) != POSITIVE){
    roots.push_back(ur);
    if (lr >lb-max_error && lr < ub){
    roots.push_back(lr);
    }
    }
    } else {
    if (lr > lb-max_error && lr <ub){
    if (!CLEAN || sign(begin[2]) != NEGATIVE){
    roots.push_back(lr);
    }
    }
    }*/
}
}

CGAL_INLINE_FUNCTION
void compute_quadratic_roots(const double *begin, const double *end,
double lb, double ub, std::vector<double> &roots)
{
    return compute_quadratic_roots_t<false>(begin, end, lb, ub, roots);
}

CGAL_INLINE_FUNCTION
void compute_quadratic_cleaned_roots(const double *begin, const double *end,
double lb, double ub, std::vector<double> &roots)
{
    return compute_quadratic_roots_t<true>(begin, end, lb, ub, roots);
}

namespace {
template <bool CLEAN, class NT>
inline void compute_linear_roots_t(const NT *begin, const NT *,
					  NT lb, NT ub,
					  std::vector<NT> &roots)
{
    if (CLEAN &&  begin[1]>0 ) return;
    //NT max_error=0;
    //if (CLEAN) max_error=max_error_value;
    NT r= -CGAL::to_double(begin[0]/begin[1]);
    if ((CLEAN || r > lb) && r < ub) {
        roots.push_back(r);
    }
}
}

CGAL_INLINE_FUNCTION
void compute_linear_roots(const double *begin, const double *end,
double lb, double ub, std::vector<double> &roots)
{
    return compute_linear_roots_t<false>(begin, end, lb, ub, roots);
}

CGAL_INLINE_FUNCTION
void compute_linear_cleaned_roots(const double *begin, const double *end,
				  double lb, double ub, std::vector<double> &roots)
{
    return compute_linear_roots_t<true>(begin, end, lb, ub, roots);
}


namespace {
template <class NT>
 inline void filter_roots_t(const NT *begin, const NT *end,
				  NT lb, NT ub, NT last_root, std::vector<NT> &roots)
{
// if we are not close to the current time, then we are fine
    if (roots.empty()) return;
    //if (roots.back() > lb+ .0005) return;

    //double eps= .0005;
    /*double last_root=-std::numeric_limits<double>::infinity();

    while (roots.back() < lb) {
      last_root= roots.back();
      roots.pop_back();
      }*/
    //if (roots.back() > lb+eps) return;

    //typedef CGAL_POLYNOMIAL_NS::Polynomial<NT> Fn;

    typedef CGAL_POLYNOMIAL_NS::Interval_polynomial IFn;
    typedef CGAL_POLYNOMIAL_NS::internal::Derivative<IFn> Diff;
    typedef typename IFn::NT INT;
    

    /*bool popped=false;*/
    // if the last valid root is closer than last, consider it as doubtful instead
    if (lb-last_root > roots.back()-lb) {
      last_root= roots.back();
      roots.pop_back();
      /*popped=true;*/
    } /*else {
      last_root=lb;
      }*/

    INT vi;
    if (last_root== -std::numeric_limits<double>::infinity()){
      if ((end-begin)%2==1) {
	vi= std::numeric_limits<double>::infinity();
      } else {
	vi = -*(end-1);
      }
    } else {
      IFn fi(begin, end);
      if (roots.empty()) {
	Interval_arithmetic_guard guard;
        if (ub== std::numeric_limits<double>::infinity()) {
          vi = 10*lb + 1000;
        } else {
          vi = fi((INT(lb)+INT(ub))/2.0);
        }
      } else {
	Interval_arithmetic_guard guard;
	vi = fi((INT(last_root)+INT(roots.back()))/2.0);
      }
    }
    
    if (vi.inf() > 0) {
      return;
    } else if (vi.sup() < 0){
      roots.push_back(last_root);
    
      /*if (!popped) {
	IFn f(begin, end);
	std::cout << "Adding last due to sign of " << vi << std::endl;
	std::cout << "last " << last_root << " lb " << lb << " poly " << f << std::endl;
	}*/
      return;
    }
    Interval_arithmetic_guard guard;
    Diff dx;
    IFn f(begin, end);
    IFn d= dx(f);

    INT dv= d(roots.back());
    // switch
    //while (sign(d(roots.back().representation()))== ZERO) d= dx_(d);
    while (dv.inf() <= 0 && dv.sup() >= 0) {
      d= dx(d);
      dv= d(roots.back());
    }
    // switch
    //if (sign(d(roots.back().representation()))==POSITIVE){
    if (dv.sup() < 0) {
        roots.push_back(last_root);
	/*if (!popped) {
	  IFn f(begin, end);
	  std::cout << "Adding last due to deriv of " << vi << std::endl;
	  std::cout << "last " << last_root << " lb " << lb << " poly " << f << std::endl;
	  }*/
    }
}
}

CGAL_INLINE_FUNCTION
void filter_solver_roots(const double *begin, const double *end,
			 double lb, double ub, double last,
			 std::vector<double> &roots)
{
  filter_roots_t(begin, end, lb, ub, last, roots);
}


/*void polynomial_compute_roots(const double *begin, const double *end,  double lb,
                  double ub, std::vector<double> &roots){
#ifdef POLYNOMIAL_USE_GSL
    gsl_polynomial_compute_roots(begin, end, lb, ub, roots);
#else
    jama_polynomial_compute_roots(begin, end, lb, ub, roots);
#endif
}

void polynomial_compute_cleaned_roots(const double *begin, const double *end,  double lb,
                      double ub, std::vector<double> &roots){
#ifdef POLYNOMIAL_USE_GSL
gsl_polynomial_compute_cleaned_roots(begin, end, lb, ub, roots);
#else
jama_polynomial_compute_cleaned_roots(begin, end, lb, ub, roots);
#endif
}*/

CGAL_INLINE_FUNCTION
double evaluate_polynomial(const double *b, const double *e, double t)
{
#ifdef POLYNOMIAL_USE_GSL
    return gsl_evaluate_polynomial(b, e, t);
#else
    if (b==e) return 0.0;

    const double *rit=e-1;
    double result = *rit;
    --rit;
    for (; rit != b-1; --rit) {
        result *= t;
        result += (*rit);
    }
    return result;
#endif
}


} } } //namespace CGAL::POLYNOMIAL::internal
