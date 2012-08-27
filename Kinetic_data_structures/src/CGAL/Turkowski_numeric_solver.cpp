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

/*#ifdef _MSC_VER
#pragma warning(disable:1572)
#endif*/


//#include "numeric_solvers_support.C"

// Taken from http://www.worldserver.com/turk/opensource/

/* Copyright (C) 1978-1999 Ken Turkowski. <turk_at_computer.org>
 *
 * All rights reserved.
 *
 * Warranty Information
 *  Even though I have reviewed this software, I make no warranty
 *  or representation, either express or implied, with respect to this
 *  software, its quality, accuracy, merchantability, or fitness for a
 *  particular purpose.  As a result, this software is provided "as is,"
 *  and you, its user, are assuming the entire risk as to its quality
 *  and accuracy.
 *
 * This code may be used and freely distributed as long as it includes
 * this copyright notice and the above warranty information.
 */

#include <cmath>
#include <iomanip>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

# define FLOAT double

static const unsigned int MAXN= 55;
//# define PARAMFLOAT double_t

/*******************************************************************************
 * FindCubicRoots
 *
 *	Solve:
 *		coeff[3] * x^3 + coeff[2] * x^2 + coeff[1] * x + coeff[0] = 0
 *
 *	returns:
 *		3 - 3 real roots
 *		1 - 1 real root (2 complex conjugate)
 *******************************************************************************/

static long
FindCubicRoots(const FLOAT coeff[4], FLOAT x[3])
{
  FLOAT a1 = coeff[2] / coeff[3];
  FLOAT a2 = coeff[1] / coeff[3];
  FLOAT a3 = coeff[0] / coeff[3];

  double Q = (a1 * a1 - 3 * a2) / 9;
  double R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
  double Qcubed = Q * Q * Q;
  double d = Qcubed - R * R;

  /* Three real roots */
  if (d >= 0) {
    double theta = std::acos(R / std::sqrt(Qcubed));
    double sqrtQ = std::sqrt(Q);
    x[0] = -2 * sqrtQ * std::cos( theta           / 3) - a1 / 3;
    x[1] = -2 * sqrtQ * std::cos((theta + 2 * CGAL_PI) / 3) - a1 / 3;
    x[2] = -2 * sqrtQ * std::cos((theta + 4 * CGAL_PI) / 3) - a1 / 3;
    return (3);
  }

  /* One real root */
  else {
    double e = std::pow(std::sqrt(-d) + ::CGAL::abs(R), 1. / 3.);
    if (R > 0)
      e = -e;
    x[0] = (e + Q / e) - a1 / 3.;
    return (1);
  }
}


/*******************************************************************************
 * FindPolynomialRoots
 *
 * The Bairstow and Newton correction formulae are used for a simultaneous
 * linear and quadratic iterated synthetic division.  The coefficients of
 * a polynomial of degree n are given as a[i] (i=0,i,..., n) where a[0] is
 * the constant term.  The coefficients are scaled by dividing them by
 * their geometric mean.  The Bairstow or Newton iteration method will
 * nearly always converge to the number of figures carried, fig, either to
 * root values or to their reciprocals.  If the simultaneous Newton and
 * Bairstow iteration fails to converge on root values or their
 * reciprocals in maxiter iterations, the convergence requirement will be
 * successively reduced by one decimal figure.  This program anticipates
 * and protects against loss of significance in the quadratic synthetic
 * division.  (Refer to "On Programming the Numerical Solution of
 * Polynomial Equations," by K. W. Ellenberger, Commun. ACM 3 (Dec. 1960),
 * 644-647.)  The real and imaginary part of each root is stated as u[i]
 * and v[i], respectively, together with the corresponding constant,
 * conv[i], used in the convergence test.  This program has been used
 * successfully for over a year on the Bendix G15-D (Intercard System) and
 * has recently been coded for the IBM 709 (Fortran system).
 *
 * ACM algorithm #30 - Numerical Solution of the Polynomial Equation
 * K. W. Ellenberger
 * Missle Division, North American Aviation, Downey, California
 * Converted to C, modified, optimized, and structured by
 * Ken Turkowski
 * CADLINC, Inc., Palo Alto, California
 *******************************************************************************/

static void
FindPolynomialRoots(
		    const FLOAT     *a,                               /* Coefficients */
		    FLOAT           *u,                               /* Real component of each root */
		    FLOAT           *v,                               /* Imaginary component of each root */
		    FLOAT           *conv,                            /* Convergence constant associated with each root */
		    register long   n,                                /* Degree of polynomial (order-1) */
		    long            maxiter,                          /* Maximum number of iterations */
		    long            fig                               /* The number of decimal figures to be computed */
		    )
{
  int number_of_ITERATE=0;
  int number_of_INIT=0;
  CGAL_precondition(static_cast<unsigned int>(fig) < MAXN);
  int i;
  register int j;
  FLOAT h[MAXN + 3], b[MAXN + 3], c[MAXN + 3], d[MAXN + 3], e[MAXN + 3];
  /* [-2 : n] */
  FLOAT K, ps, qs, pt, qt, s, rev, r= std::numeric_limits<double>::infinity();
  int t;
  FLOAT p=std::numeric_limits<double>::infinity(), q=std::numeric_limits<double>::infinity();

  /* Zero elements with negative indices */
  b[2 + -1] = b[2 + -2] =
    c[2 + -1] = c[2 + -2] =
    d[2 + -1] = d[2 + -2] =
    e[2 + -1] = e[2 + -2] =
    h[2 + -1] = h[2 + -2] = 0.0;

  /* Copy polynomial coefficients to working storage */
  for (j = n; j >= 0; j--)
    h[2 + j] = *a++;                          /* Note reversal of coefficients */

  t = 1;
  K = std::pow(10.0, (double)(fig));            /* Relative accuracy */

  for (; h[2 + n] == 0.0; n--) {                /* Look for zero high-order coeff. */
    *u++ = 0.0;
    *v++ = 0.0;
    *conv++ = K;
  }

 INIT:
  ++number_of_INIT;
  if (number_of_INIT > 1000) {
    std::cerr << "Too many INITs" << std::flush;
    return;
  }

  if (n == 0)
    return;

  ps = qs = pt = qt = s = 0.0;
  rev = 1.0;
  K = std::pow(10.0, (double)(fig));

  if (n == 1) {
    r = -h[2 + 1] / h[2 + 0];
    goto LINEAR;
  }

  for (j = n; j >= 0; j--)                      /* Find geometric mean of coeff's */
    if (h[2 + j] != 0.0)
      s += std::log( ::CGAL::abs(h[2 + j]));
  s = std::exp(s / (n + 1));

  for (j = n; j >= 0; j--)                      /* Normalize coeff's by mean */
    h[2 + j] /= s;

  if ( ::CGAL::abs(h[2 + 1] / h[2 + 0]) < ::CGAL::abs(h[2 + n - 1] / h[2 + n])) {
  REVERSE:
    t = -t;
    for (j = (n - 1) / 2; j >= 0; j--) {
      s = h[2 + j];
      h[2 + j] = h[2 + n - j];
      h[2 + n - j] = s;
    }
  }
  if (qs != 0.0) {
    p = ps;
    q = qs;
  }
  else {
    if (h[2 + n - 2] == 0.0) {
      q = 1.0;
      p = -2.0;
    }
    else {
      q = h[2 + n] / h[2 + n - 2];
      p = (h[2 + n - 1] - q * h[2 + n - 3]) / h[2 + n - 2];
    }
    if (n == 2)
      goto QADRTIC;
    r = 0.0;
  }
 ITERATE:
  ++number_of_ITERATE;
  if (number_of_ITERATE > 1000) {
    std::cerr << "Too many ITERATEs" << std::flush;
    return;
  }
  for (i = maxiter; i > 0; i--) {

    for (j = 0; j <= n; j++) {                /* BAIRSTOW */
      b[2 + j] = h[2 + j] - p * b[2 + j - 1] - q * b[2 + j - 2];
      c[2 + j] = b[2 + j] - p * c[2 + j - 1] - q * c[2 + j - 2];
    }
    if ((h[2 + n - 1] != 0.0) && (b[2 + n - 1] != 0.0)) {
      if ( ::CGAL::abs(h[2 + n - 1] / b[2 + n - 1]) >= K) {
	b[2 + n] = h[2 + n] - q * b[2 + n - 2];
      }
      if (b[2 + n] == 0.0)
	goto QADRTIC;
      if (K < ::CGAL::abs(h[2 + n] / b[2 + n]))
	goto QADRTIC;
    }

    for (j = 0; j <= n; j++) {                /* NEWTON */
      /* Calculate polynomial at r */
      d[2 + j] = h[2 + j] + r * d[2 + j - 1];
      /* Calculate derivative at r */
      e[2 + j] = d[2 + j] + r * e[2 + j - 1];
    }
    if (d[2 + n] == 0.0)
      goto LINEAR;
    if (K < ::CGAL::abs(h[2 + n] / d[2 + n]))
      goto LINEAR;

    c[2 + n - 1] = -p * c[2 + n - 2] - q * c[2 + n - 3];
    s = c[2 + n - 2] * c[2 + n - 2] - c[2 + n - 1] * c[2 + n - 3];
    if (s == 0.0) {
      p -= 2.0;
      q *= (q + 1.0);
    }
    else {
      p += (b[2 + n - 1] * c[2 + n - 2] - b[2 + n] * c[2 + n - 3]) / s;
      q += (-b[2 + n - 1] * c[2 + n - 1] + b[2 + n] * c[2 + n - 2]) / s;
    }
    if (e[2 + n - 1] == 0.0)
      r -= 1.0;                             /* Minimum step */
    else
      r -= d[2 + n] / e[2 + n - 1];         /* Newton's iteration */
  }
  ps = pt;
  qs = qt;
  pt = p;
  qt = q;
  if (rev < 0.0)
    K /= 10.0;
  rev = -rev;
  goto REVERSE;

 LINEAR:
  if (t < 0)
    r = 1.0 / r;
  n--;
  *u++ = r;
  *v++ = 0.0;
  *conv++ = K;

  for (j = n; j >= 0; j--) {                    /* Polynomial deflation by lin-nomial */
    if ((d[2 + j] != 0.0) && ( ::CGAL::abs(h[2 + j] / d[2 + j]) < K))
      h[2 + j] = d[2 + j];
    else
      h[2 + j] = 0.0;
  }

  if (n == 0)
    return;
  goto ITERATE;

 QADRTIC:
  if (t < 0) {
    p /= q;
    q = 1.0 / q;
  }
  n -= 2;

  if (0.0 < (q - (p * p / 4.0))) {              /* Two complex roots */
    *(u + 1) = *u = -p / 2.0;
    u += 2;
    s = sqrt(q - (p * p / 4.0));
    *v++ = s;
    *v++ = -s;
  }                                             /* Two real roots */
  else {
    s = std::sqrt(((p * p / 4.0)) - q);
    if (p < 0.0)
      *u++ = -p / 2.0 + s;
    else
      *u++ = -p / 2.0 - s;
    *u = q / u[-1];
    ++u;                                      // moved from lhs of before
    *v++ = 0.0;
    *v++ = 0.0;
  }
  *conv++ = K;
  *conv++ = K;

  for (j = n; j >= 0; j--) {                    /* Polynomial deflation by quadratic */
    if ((b[2 + j] != 0.0) && ( ::CGAL::abs(h[2 + j] / b[2 + j]) < K))
      h[2 + j] = b[2 + j];
    else
      h[2 + j] = 0.0;
  }
  goto INIT;
}


#undef MAXN

template <bool CLEAN>
static void Turkowski_polynomial_compute_roots_t(const double *begin, 
						 const double *end,
						 double lb, double ub,
						 std::vector<double> &roots)
{
  std::size_t numc= end-begin;
  double rp[MAXN];
  double cp[MAXN];
  double cc[MAXN];

  /*for (unsigned int i=0; i< numc; ++i){
    rp[i]= std::numeric_limits<double>::infinity();
    }*/

  for (unsigned int i=0; i< MAXN; ++i){
    cc[i]=std::numeric_limits<double>::infinity();
    rp[i]=std::numeric_limits<double>::infinity();
    cc[i]=std::numeric_limits<double>::infinity();
  }

  FindPolynomialRoots(begin, rp, cp, cc, static_cast<long>(numc)-1, 10*static_cast<long>(numc), 40);

  /*if (CLEAN) {
    lb-= .000005;
    }*/
  double last= -std::numeric_limits<double>::infinity();
  for (std::size_t i=0; i< numc-1; ++i) {
    /* std::cout << "Trying " <<  rp[i] << "+" << std::setprecision(10) <<  cp[i] << "i "
       << cc[i] << "\n";*/
    if (cc[i] > 10000 && root_is_good(rp[i], cp[i], lb, ub)) {
      roots.push_back(rp[i]);
      /*std::cout << "Good was " <<  rp[i] << "+" <<std::setprecision(10) <<  cp[i] << "i "
	<< cc[i] << "\n";*/
    } else if (CLEAN && root_is_good(rp[i], cp[i], last, ub) && last < rp[i] && rp[i] <= lb) {
      last= rp[i];
    }
    /*else {
      std::cout << "Rejected " << rp[i] << "+" << cp[i] << "i\n";
      }*/
  }
  std::sort(roots.begin(), roots.end(), std::greater<double>());
  
  if (CLEAN) filter_solver_roots(begin, end, lb, ub, last, roots);
}


void Turkowski_polynomial_compute_roots(const double *begin, const double *end, 
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
  case 3:
    {
      double rd[3];
      int numr= FindCubicRoots(begin, rd);
      for (int i=numr-1; i>=0; --i) {
	if (rd[i] >= lb && rd[i] < ub) roots.push_back(rd[i]);
      }
      std::sort(roots.begin(), roots.end(), std::greater<double>());
      break;
    }
  default:
    Turkowski_polynomial_compute_roots_t<false>(begin, end, lb, ub, roots);

  }

}


void Turkowski_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
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
  case 3:
    {
      double rd[3];
      int numr= FindCubicRoots(begin, rd);
      double last=-std::numeric_limits<double>::infinity();
      for (int i=numr-1; i>=0; --i) {
	if (rd[i]< ub && rd[i] >= lb) roots.push_back(rd[i]);
	if (rd[i] < lb && rd[i] > last){
	  last=rd[i];
	}
      }
      std::sort(roots.begin(), roots.end(), std::greater<double>());
      filter_solver_roots(begin, end, lb, ub, last, roots);
      break;
    }
  default:
    Turkowski_polynomial_compute_roots_t<true>(begin, end, lb, ub, roots);
  }
}


} } } //namespace CGAL::POLYNOMIAL::internal
