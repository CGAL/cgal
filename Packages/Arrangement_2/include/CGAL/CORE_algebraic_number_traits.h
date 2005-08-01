// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_CORE_ALGEBRAIC_NUMBER_TRAITS_2_H
#define CGAL_CORE_ALGEBRAIC_NUMBER_TRAITS_2_H

/*! \file
 * The number-type traits for CORE algebraic numbers.
 */

#include <CORE/BigInt.h>
#include <CORE/BigRat.h>
#include <CGAL/CORE_Expr.h>

CGAL_BEGIN_NAMESPACE

/*!
 * \class A traits class for CORE's algebraic number types.
 */
class CORE_algebraic_number_traits
{
public:

  typedef CORE::BigInt                    Integer;
  typedef CORE::BigRat                    Rational;
  typedef CORE::Polynomial<Integer>       Polynomial;
  typedef CORE::Expr                      Algebraic;

  /*!
   * Convert an integer to an algebraic number.
   * \param z An integer.
   * \return An algebraic number equivalent to z.
   */
  Algebraic convert (const Integer& z) const
  {
    return (Algebraic (z));
  }

  /*!
   * Convert a rational number to an algebraic number.
   * \param q A rational number.
   * \return An algebraic number equivalent to q.
   */
  Algebraic convert (const Rational& q) const
  {
    return (Algebraic (q));
  }

  /*!
   * Convert a sequence of rational coefficients to an equivalent sequence
   * of integer coefficients. If the input coefficients are q(1), ..., q(k),
   * where q(i) = n(i)/d(i) then the output coefficients will be of the
   * form:
   *               n(i) * lcm {d(1), ... , d(k)} 
   *       a(i) = -------------------------------
   *               d(i) * gcd {n(1), ... , n(k)}
   *
   * \param q_begin The begin iterator of the rational sequence.
   * \param q_end The past-the-end iterator of the rational sequence.
   * \param zoi An output iterator for the integer coefficients.
   * \return A past-the-end iterator for the output container.
   * \pre The value type of q_begin and q_end is Rational,
   *      and the value type of zoi is Integer.
   */
  template <class InputIterator, class OutputIterator>
  OutputIterator convert_coefficients (InputIterator q_begin,
                                       InputIterator q_end,
                                       OutputIterator zoi) const
  {
    // Compute the least common multiplicand (LCM) of the denominators,
    // denoted L, and the greatest common divisor (GCD) of the numerators,
    // denoted G.
    InputIterator  q_iter = q_begin;
    Integer        denom_lcm, temp_lcm;
    Integer        numer_gcd, temp_gcd;
    Integer        numer, denom;

    denom_lcm = CORE::denominator (*q_iter);
    numer_gcd = CORE::numerator (*q_iter);

    ++q_iter;
    while (q_iter != q_end)
    {
      numer = CORE::numerator (*q_iter);

      if (CGAL::sign (numer) != ZERO)
      {
        denom = CORE::denominator (*q_iter);

        temp_lcm = denom_lcm;
        temp_gcd = numer_gcd;

        denom_lcm *= denom;
        denom_lcm /= CORE::gcd (temp_lcm, denom);

        numer_gcd = CORE::gcd (temp_gcd, numer);
      }

      ++q_iter;
    }

    // Generate the output coefficients (n(i)*L) / (d(i)*G).
    for (q_iter = q_begin; q_iter != q_end; ++q_iter)
    {
      *zoi = (CORE::numerator (*q_iter) * denom_lcm) /
             (numer_gcd * CORE::denominator (*q_iter));

      ++zoi;
    }

    return (zoi);
  }

  /*!
   * Compute the square root of an algebraic number.
   * \param x The number.
   * \return The sqaure root of x.
   * \pre x is non-negative.
   */
  Algebraic sqrt (const Algebraic& x) const
  {
    return (CGAL::sqrt (x));
  }

  /*!
   * Compute the roots of a quadratic equations with integer coefficients:
   *   a*x^2+ b*x + c = 0
   * \param oi An output iterator for the real-valued solutions of the
   *           quadratic equation.
   * \return A past-the-end iterator for the output container.
   * \pre The value type of oi is Algebraic.
   */
  template <class NT, class OutputIterator>
  OutputIterator solve_quadratic_equation (const NT& a,
                                           const NT& b,
                                           const NT& c,
                                           OutputIterator oi) const
  {
    // Check if this is really a linear equation.
    const Sign     sign_a = CGAL::sign (a);

    if (sign_a == ZERO)
    {
      // Solve a linear equation.
      if (CGAL::sign(b) != ZERO)
      {
        *oi = -Algebraic (c) / Algebraic (b);
        ++oi;
      }

      return (oi);
    }

    // Act according to the discriminant.
    const NT       disc = b*b - 4*a*c;
    const Sign     sign_disc = CGAL::sign (disc);


    if (sign_disc == ZERO)
    {
      // We have one real root with mutliplicity 2.
      *oi = -Algebraic (b) / Algebraic (2*a);
      ++oi;
    }
    else if (sign_disc == POSITIVE)
    {
      // We have two distinct real roots. We return them in ascending order.
      const Algebraic      sqrt_disc = CGAL::sqrt (Algebraic (disc));
      const Algebraic      alg_b = b;
      const Algebraic      alg_2a = 2*a;

      if (sign_a == POSITIVE)
      {
        *oi = -(sqrt_disc + alg_b) / alg_2a;
        ++oi;
        *oi = (sqrt_disc - alg_b) / alg_2a;
        ++oi;
      }
      else
      {
        *oi = (sqrt_disc - alg_b) / alg_2a;
        ++oi;
        *oi = -(sqrt_disc + alg_b) / alg_2a;
        ++oi;
      }
    }

    return (oi);
  }

  /*!
   * Construct a polynomial with integer coefficients.
   * \param coeffs The coefficients of the input polynomial.
   * \param degree The degree of the input polynomial.
   * \return The polynomial.
   */
  Polynomial construct_polynomial (const Integer *coeffs,
				   unsigned int degree) const
  {
    Polynomial   poly = Polynomial (degree, const_cast<Integer*> (coeffs));
    poly.contract();

    return (poly);
  }

  /*!
   * Construct a polynomial with integer coefficients given rational
   * coefficients.
   * \param coeffs The coefficients of the input polynomial.
   * \param degree The degree of the input polynomial.
   * \param poly Output: The resulting polynomial with integer coefficients.
   * \param poly_denom Output: The denominator for the polynomial.
   * \return Whether this polynomial is non-zero (false if the polynomial is 
   *         zero).
   */
  bool construct_polynomial (const Rational *coeffs,
			     unsigned int degree,
			     Polynomial& poly, Integer& poly_denom) const
  {
    // Compute the true degree for the polynomial.
    while (CGAL::sign (coeffs[degree]) == ZERO)
    {
      if (degree == 0)
	return (false);
      degree--;
    }

    // Compute the least common multiplicand (LCM) of the denominators,
    // denoted L, and the greatest common divisor (GCD) of the numerators,
    // denoted G.
    int            index = degree;
    Integer        denom_lcm, temp_lcm;
    Integer        numer_gcd, temp_gcd;
    Integer        numer, denom;

    denom_lcm = CORE::denominator (coeffs[index]);
    numer_gcd = CORE::numerator (coeffs[index]);

    index--;
    while (index >= 0)
    {
      numer = CORE::numerator (coeffs[index]);

      if (CGAL::sign (numer) != ZERO)
      {
        denom = CORE::denominator (coeffs[index]);

        temp_lcm = denom_lcm;
        temp_gcd = numer_gcd;

        denom_lcm *= denom;
        denom_lcm /= CORE::gcd (temp_lcm, denom);

        numer_gcd = CORE::gcd (temp_gcd, numer);
      }

      index--;
    }

    // Generate the output coefficients (n(i)*L) / (d(i)*G).
    Integer                *z_coeffs = new Integer [degree + 1];

    for (index = 0; index <= static_cast<int>(degree); index++)
    {
      z_coeffs[index] = (CORE::numerator (coeffs[index]) * denom_lcm) /
                        (numer_gcd * CORE::denominator (coeffs[index]));
    }

    // Set the output.
    poly = Polynomial (degree, z_coeffs);
    poly_denom = denom_lcm;

    delete[] z_coeffs;
    return (true);
  }

  /*!
   * Scale a polynomial by a given factor.
   * \param poly Input: The polynomial.
   *             Output: The scaled polynomial.
   * \param factor The scaling factor.
   * \return The scaled polynomial.
   */
  void scale (Polynomial& poly, const Integer& factor) const
  {
    poly.mulScalar (factor);
    return;
  }

  /*!
   * Compute the degree of a polynomial.
   */
  int degree (Polynomial& poly) const
  {
    return (poly.getTrueDegree());
  }

  /*!
   * Evaluate a polynomial at a given x-value.
   * \param x The value to evaluate at.
   * \return The value of the polynomial at x.
   */
  template <class NT>
  NT evaluate_at (const Polynomial& poly,
		  NT& x) const
  {
    return (poly.eval (x));
  }

  /*!
   * Compute the derivative of the given polynomial.
   * \param poly The polynomial p(x).
   * \return The derivative p'(x).
   */
  Polynomial derive (const Polynomial& poly) const
  {
    return (differentiate (poly));
  }

  /*!
   * Compute the real-valued roots of a polynomial with integer coefficients,
   * sorted in ascending order.
   * \param poly The input polynomial.
   * \param oi An output iterator for the real-valued root of the polynomial.
   * \return A past-the-end iterator for the output container.
   * \pre The value type of oi is Algebraic.
   */
  template <class OutputIterator>
  OutputIterator compute_polynomial_roots (Polynomial& poly,
					   OutputIterator oi) const
  {
    // Get the real degree of the polynomial.
    int        degree = poly.getTrueDegree();

    if (degree <= 0)
      return (oi);

    // Check if we really have a simple quadratic equation.
    if (degree <= 2)
    {
      return (solve_quadratic_equation ((degree == 2 ? poly.getCoeff(2) : 0), 
					poly.getCoeff(1),
					poly.getCoeff(0),
					oi));
    }

    // Compute the real-valued roots of the polynomial.
    CORE::Sturm<Integer>       sturm (poly);
    const unsigned int         n_roots = sturm.numberOfRoots();
    unsigned int               i;
    
    for (i = 1; i <= n_roots; i++)
    {
      // Get the i'th real-valued root.
      *oi = rootOf(poly, i);
      ++oi;
    }

    return (oi);
  }
};

CGAL_END_NAMESPACE

#endif
