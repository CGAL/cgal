// Copyright (c) 1999  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_CONIC_ARC_2_EQ_CORE_H
#define CGAL_CONIC_ARC_2_EQ_CORE_H

CGAL_BEGIN_NAMESPACE

/*!
 * Solve the quadratic equation: a*x^2 + b*x + c = 0.
 * Note that the number types of the solutions and of the equation coefficients
 * need not be the same, but it should be possible to construct a SolNT number
 * from a CfNT.
 * \param a The coefficient of x^2.
 * \param b The coefficient of x.
 * \param c The free coefficient.
 * \param roots The solutions of the equation.
 *              This area must be allocated to the size of 2.
 * \return The number of distinct solutions.
 */
template <class CfNT, class SolNT>
int solve_quadratic_eq (const CfNT& a, const CfNT& b, const CfNT& c,
			SolNT* roots)
{
  // Check if this is really a linear equation.
  const CfNT _zero = 0;

  if (a == _zero)
  {
    // Solve a linear equation.
    if (b != _zero)
    {
      roots[0] = -SolNT(c) / SolNT(b);
      return (1);
    }
    else
    {
      return (0);
    }
  }

  // Act according to the discriminant.
  const CfNT        _two = 2;
  const CfNT        _four = 4;
  const CfNT        disc = b*b - _four*a*c;
  Comparison_result res = CGAL_NTS compare(disc, _zero);

  if (res == SMALLER)
  {
    // No real roots.
    return (0);
  }
  else if (disc == EQUAL)
  {
    // One real root with mutliplicity 2.
    roots[0] = -SolNT(b) / SolNT(_two*a);
    return (1);
  }
  else
  {
    // Two real roots.
    const SolNT      sqrt_disc = CGAL::sqrt(SolNT(disc));
 
    roots[0] = (sqrt_disc - SolNT(b)) / SolNT(_two*a);
    roots[1] = -(sqrt_disc + SolNT(b)) / SolNT(_two*a);
    return (2);
  }
}

/*!
 * Solve a quartic equation: a*x^4 + b*x^3 + c*x^2 + d*x + e = 0.
 * Note that the number types of the solutions and of the equation coefficients
 * need not be the same, but it should be possible to construct a SolNT number
 * from a CfNT.
 * \param a The coefficient of x^4.
 * \param b The coefficient of x^3.
 * \param c The coefficient of x^2.
 * \param d The coefficient of x.
 * \param e The free coefficient.
 * \param roots The solutions of the equation.
 *              This area must be allocated to the size of 4.
 * \return The number of distinct solutions.
 */
template <class CfNT, class SolNT>
int solve_quartic_eq (const CfNT& a, const CfNT& b, const CfNT& c, 
		      const CfNT& d, const CfNT& e,
		      SolNT* roots)
{
  // First check whether we have 0 as a multiple root.
  static const CfNT _zero = 0;

  if (e == _zero && d == _zero)
  {
    if (c == _zero)
    {
      if (b == _zero)
      {
	// Unless the equation is trivial, we have the root 0,
	// with multiplicity of 4.
	if (a == _zero)
	  return (0);
	  
	roots[0] = SolNT(_zero);
	return (1);
      }
      else
      {
	// Add the solution 0 (with multiplicity 3), and add the solution
	// to a*x + b = 0.
	if (a == _zero)
	  return (0);
	
	roots[0] = SolNT(_zero);
	roots[1] = -SolNT(b) / SolNT(a);
	return (2);
      }
    }
    else
    {
      // Add the solution 0 (with multiplicity 2),
      // and solve a*x^2 + b*x + c*x = 0.
      roots[0] = SolNT(_zero);

      return (1 + solve_quadratic_eq<CfNT,SolNT> (a, b, c,
						  roots + 1));
    }
  }

  // In case we have an equation of a lower degree:
  if (a == _zero && b == _zero)
  {
    // We have to solve c*x^2 + d*x + e = 0.
    return (solve_quadratic_eq<CfNT, SolNT> (c, d, e,
					     roots));
  }

  // In case we have an equation of the form:
  //  a*x^4 + c*x^2 + e = 0
  //
  // Then by substituting y = x^2, we obtain a quadratic equation:
  //  a*y^2 + c*y + e = 0
  //
  if (b == _zero && d == _zero)
  {
    // Solve the equation for y = x^2.
    SolNT   sq_roots[2];
    int     n_sq;

    n_sq = solve_quadratic_eq<CfNT, SolNT> (a, c, e,
					    sq_roots);

    // Convert to roots of the original equation: only positive roots for y
    // are relevant of course.
    int  n = 0;
    int  j;

    for (j = 0; j < n_sq; j++)
    {
      if (sq_roots[j] > SolNT(0))
      {
	roots[n] = CGAL::sqrt(sq_roots[j]);
	roots[n+1] = -roots[n];
	n += 2;
      }
    }

    return (n);
  }

  // Set a quadratic polynomial p(x) and compute its real-valued roots.
  std::vector<CfNT>       coeffs(5);

  coeffs[4] = a;
  coeffs[3] = b;
  coeffs[2] = c;
  coeffs[1] = d;
  coeffs[0] = e;

  CORE::Polynomial<CfNT> p (coeffs);
  CORE::Sturm<CfNT>      sturm (p);
  const int              n_roots = sturm.numberOfRoots();
  int                    i;

  for (i = 1; i <= n_roots; i++)
  {
    // Get the i'th real-valued root.
    roots[i - 1] = rootOf(p, i);
  }

  return (n_roots);
}

CGAL_END_NAMESPACE

#endif
