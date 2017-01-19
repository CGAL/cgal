// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_CORE_ALGEBRAIC_NUMBER_TRAITS_2_H
#define CGAL_CORE_ALGEBRAIC_NUMBER_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * The number-type traits for CORE algebraic numbers.
 */

#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_Expr.h>

namespace CGAL {

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
   * Get the numerator of a rational number.
   * \param q A rational number.
   * \return The numerator of q.
   */
  Integer numerator (const Rational& q) const
  {
    return (CORE::numerator (q));
  }

  /*!
   * Get the denominator of a rational number.
   * \param q A rational number.
   * \return The denominator of q.
   */
  Integer denominator (const Rational& q) const
  {
    return (CORE::denominator (q));
  }

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
   * Construct a rational number that lies strictly between two algebraic
   * values.
   * \param x1 The first algebraic value.
   * \param x2 The second algebraic value.
   * \pre The two values are not equal.
   * \return A rational number that lies in the open interval (x1, x2).
   */
  Rational rational_in_interval (const Algebraic& x1,
                                 const Algebraic& x2) const
  {
    CGAL_precondition (x1 != x2);

    const Integer  one (1);
    Algebraic      scaled_x1 = x1;
    Integer        ix1;
    Algebraic      scaled_x2 = x2;
    Integer        ix2;
    Integer        denom = 1;

    while (true)
    {
      // Check if x1*2^k and x2*2^k are separated by at least 2.
      ix1 = scaled_x1.BigIntValue();
      ix2 = scaled_x2.BigIntValue();

      if (CORE::abs (ix2 - ix1) > one)
        break;

      // Scale the values by a factor of 2.
      scaled_x1 *= 2;
      scaled_x2 *= 2;
      denom *= 2;
    }

    // If we reached a separation, we construct a rational that
    // approximates (x1 + x2)/2.
    return (Rational (ix1 + ix2, 2*denom));
  }

  /*!
   * Get a range of double-precision floats that contain the given algebraic
   * number.
   * \param x The given number.
   * \return A pair <x_lo, x_hi> that contain x.
   */
  std::pair<double, double> double_interval (const Algebraic& x) const
  {
    double      x_lo, x_hi;

    x.doubleInterval (x_lo, x_hi);
    return (std::make_pair (x_lo, x_hi));
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
      {
        poly_denom = 1;
	return (false);
      }
      degree--;
    }

    // Compute the least common multiplicand (LCM) of the denominators,
    // denoted L.
    int            index = degree;
    Integer        denom_lcm, temp_lcm;
    Integer        denom;

    denom_lcm = CORE::denominator (coeffs[index]);

    index--;
    while (index >= 0)
    {
      if (CGAL::sign (CORE::numerator (coeffs[index])) != ZERO)
      {
        denom = CORE::denominator (coeffs[index]);

        temp_lcm = denom_lcm;

        denom_lcm *= denom;
        denom_lcm /= CORE::gcd (temp_lcm, denom);
      }

      index--;
    }

    // Generate the output coefficients (n(i)*L) / d(i).
    Integer                *z_coeffs = new Integer [degree + 1];

    for (index = 0; index <= static_cast<int>(degree); index++)
    {
      z_coeffs[index] = (CORE::numerator (coeffs[index]) * denom_lcm) /
                        CORE::denominator (coeffs[index]);
    }

    // Set the output.
    poly = Polynomial (degree, z_coeffs);
    poly_denom = denom_lcm;

    delete[] z_coeffs;
    return (true);
  }

  /*!
   * Construct two polynomials with integer coefficients such that P(x)/Q(x)
   * is a rational function equivalent to the one represented by the two
   * given vectors of rational coefficients. It is guaranteed that the GCD
   * of P(x) and Q(x) is trivial.
   * \param p_coeffs The coefficients of the input numerator polynomial.
   * \param p_degree The degree of the input numerator polynomial.
   * \param q_coeffs The coefficients of the input denominator polynomial.
   * \param q_degree The degree of the input denominator polynomial.
   * \param p_poly Output: The resulting numerator polynomial with integer
   *                       coefficients.
   * \param q_poly Output: The resulting denominator polynomial with integer
   *                       coefficients.
   * \return (true) on success; (false) if the denominator is 0.
   */
  bool construct_polynomials (const Rational *p_coeffs,
                              unsigned int p_degree,
                              const Rational *q_coeffs,
                              unsigned int q_degree,
                              Polynomial& p_poly, Polynomial& q_poly) const
  {
    typedef CORE::Polynomial<Rational>            Rat_polynomial;

    // Construct two rational polynomials.
    Rat_polynomial   P = Rat_polynomial (p_degree,
                                         const_cast<Rational*> (p_coeffs));
    Rat_polynomial   Q = Rat_polynomial (q_degree,
                                         const_cast<Rational*> (q_coeffs));
    
    P.contract();
    Q.contract();

    int              p_deg = P.getTrueDegree();
    int              q_deg = Q.getTrueDegree();

    if (q_deg < 0)
      return (false);             // Zero denominator.

    // Check the case of a zero numerator.
    if (p_deg < 0)
    {
      Integer          coeff = 0;
      p_poly = construct_polynomial (&coeff, 0);

      coeff = 1;
      q_poly =  construct_polynomial (&coeff, 0);
      
      return (true);
    }

    // Compute the GCD of the two polynomials and normalize them.
    Rat_polynomial   g = CORE::gcd (P, Q);
    
    if (g.getTrueDegree() > 0)
    {
      P = P.pseudoRemainder (g);
      p_deg = P.getTrueDegree();
      Q = Q.pseudoRemainder (g);
      q_deg = Q.getTrueDegree();
    }

    // Construct two polynomials with integer coefficients.
    Integer        p_scale, q_scale;

    construct_polynomial (*(P.getCoeffs()), p_deg,
			  p_poly, q_scale);

    construct_polynomial (*(Q.getCoeffs()), q_deg,
			  q_poly, p_scale);

    // Scale the result polynomials.
    p_poly.mulScalar (p_scale);
    q_poly.mulScalar (q_scale);
   
    return (true);
  }

  /*!
   * Compute the degree of a polynomial.
   */
  int degree (const Polynomial& poly) const
  {
    return (poly.getTrueDegree());
  }

  /*!
   * Get the coefficient of the monomial x^i in the given polynomial.
   * \param poly The polynomial.
   * \param i The coefficient index.
   */
  Integer get_coefficient (const Polynomial& poly, unsigned int i) const
  {
    if (static_cast<int> (i) > degree(poly))
      return (0);

    return (poly.getCoeff (i));
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
   * Multiply a polynomial by some scalar coefficient.
   * \param poly The polynomial P(x).
   * \param a The scalar value.
   * \return The scalar multiplication a*P(x).
   */
  Polynomial scale (const Polynomial& poly,
                    const Integer& a) const
  {
    Polynomial   temp = poly;
    return (temp.mulScalar (a));
  }
                     
  /*!
   * Perform "long division" of two polynomials: Given A(x) and B(x) compute
   * two polynomials Q(x) and R(x) such that: A(x) = Q(x)*B(x) + R(x) and
   * R(x) has minimal degree.
   * \param polyA The first polynomial A(x).
   * \param polyB The second polynomial A(x).
   * \param rem Output: The remainder polynomial R(x).
   * \return The quontient polynomial Q(x).
   */
  Polynomial divide (const Polynomial& polyA,
                     const Polynomial& polyB,
                     Polynomial& rem) const
  {
    rem = polyA;
    return (rem.pseudoRemainder (polyB));
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
  OutputIterator compute_polynomial_roots (const Polynomial& poly,
					   OutputIterator oi) const
  {
    // Get the real degree of the polynomial.
    int            degree = poly.getTrueDegree();

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

  /*!
   * Compute the real-valued roots of a polynomial with integer coefficients,
   * within a given interval. The roots are sorted in ascending order.
   * \param poly The input polynomial.
   * \param x_min The left bound of the interval.
   * \param x_max The right bound of the interval.
   * \param oi An output iterator for the real-valued root of the polynomial.
   * \return A past-the-end iterator for the output container.
   * \pre The value type of oi is Algebraic.
   */
  template <class OutputIterator>
  OutputIterator compute_polynomial_roots (const Polynomial& poly,
                                           double x_min, double x_max,
					   OutputIterator oi) const
  {
    // Get the real degree of the polynomial.
    int            degree = poly.getTrueDegree();
    unsigned int   i;

    if (degree <= 0)
      return (oi);

    // Check if we really have a simple quadratic equation.
    if (degree <= 2)
    {
      Algebraic     alg_min (x_min), alg_max (x_max);
      Algebraic     buffer[2];
      Algebraic    *end_buffer =
        solve_quadratic_equation ((degree == 2 ? poly.getCoeff(2) : 0), 
                                  poly.getCoeff(1),
                                  poly.getCoeff(0),
                                  buffer);
      unsigned int  num_of_roots = std::distance(&buffer[0], end_buffer);

      for (i = 0; i < num_of_roots; ++i)
      {
        if (buffer[i] >= alg_min && buffer[i] <= alg_max)
        {
          *oi = buffer[i];
          ++oi;
        }
      }
      return (oi);
    }

    // Compute the real-valued roots of the polynomial in [x_min, x_max].
    CORE::Sturm<Integer>  sturm (poly);
    CORE::BFVecInterval   root_intervals;

    sturm.isolateRoots (CORE::BigFloat(x_min), CORE::BigFloat(x_max),
                        root_intervals);

    for (i = 0; i < root_intervals.size(); i++)
    {
      // Get the i'th real-valued root.
      *oi = rootOf(poly, root_intervals[i]);
      ++oi;
    }

    return (oi);
  }
};

} //namespace CGAL

#endif
