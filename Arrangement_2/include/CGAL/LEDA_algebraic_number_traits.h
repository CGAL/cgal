// Copyright (c) 2006  Barcelona.
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
//
// Author(s)     : Marcel Janer <marcel_janer@terra.es> 
//                 Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_LEDA_ALGEBRAIC_NUMBER_TRAITS_2_H
#define CGAL_LEDA_ALGEBRAIC_NUMBER_TRAITS_2_H

/*! \file
 * The number-type traits for LEDA algebraic numbers.
 */

#include <CGAL/leda_real.h>
#include <LEDA/numbers/integer.h>
#include <LEDA/numbers/rational.h> 
#include <LEDA/numbers/real.h> 
#include <LEDA/numbers/polynomial.h>
#include <LEDA/numbers/isolating.h>

#include <LEDA/core/array.h>
#include <LEDA/core/growing_array.h>

#undef LEDA_VECTOR
#define LEDA_VECTOR leda::growing_array

CGAL_BEGIN_NAMESPACE

/*!
 * \class A traits class for LEDA's algebraic number types.
 */
class LEDA_algebraic_number_traits
{
public:

  typedef leda::integer					Integer;
  typedef leda::rational				Rational;
  typedef leda::int_Polynomial			        Polynomial;
  typedef leda::real					Algebraic;

  /*!
   * Get the numerator of a rational number.
   * \param q A rational number.
   * \return The numerator of q.
   */
  Integer numerator (const Rational& q) const
  {
    return (q.numerator());
  }

  /*!
   * Get the denominator of a rational number.
   * \param q A rational number.
   * \return The denominator of q.
   */
  Integer denominator (const Rational& q) const
  {
    return (q.denominator());
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
    return (leda::small_rational_between (x1, x2));
  }

  /*!
   * Get a range of double-precision floats that contain the given algebraic
   * number.
   * \param x The given number.
   * \return A pair <x_lo, x_hi> that contain x.
   */
  std::pair<double, double> double_interval (const Algebraic& x) const
  {
    double      x_lo = x.get_lower_bound().to_double (leda::TO_N_INF);
    double      x_hi = x.get_upper_bound().to_double (leda::TO_P_INF);

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
    
    denom_lcm = (*q_iter).denominator();
    numer_gcd = (*q_iter).numerator();
    
    ++q_iter;
    while (q_iter != q_end)
    {
      numer = (*q_iter).numerator();
      
      if (leda::sign (numer) != ZERO)
      {
        denom = (*q_iter).denominator();

        temp_lcm = denom_lcm;
        temp_gcd = numer_gcd;
        
        denom_lcm *= denom;
        denom_lcm /= leda::gcd (temp_lcm, denom);
        
        numer_gcd = leda::gcd (temp_gcd, numer);
      }
      ++q_iter;
    }
    
    // Generate the output coefficients (n(i)*L) / (d(i)*G).
    for (q_iter = q_begin; q_iter != q_end; ++q_iter)
    {
      *zoi = ((*q_iter).numerator() * denom_lcm) /
        (numer_gcd * (*q_iter).denominator());
      
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
    return (leda::sqrt (x));
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
    const int     sign_a = leda::sign (a);

    if (sign_a == 0)
    {
      // Solve a linear equation.
      if (leda::sign(b) != 0)
      {
        *oi = -Algebraic (c) / Algebraic (b);
        ++oi;
      }

      return (oi);
    }

    // Act according to the discriminant.
    const NT      disc = b*b - 4*a*c;
    const int     sign_disc = leda::sign (disc);

    if (sign_disc == 0)
    {
      // We have one real root with multiplicity 2.
      *oi = -Algebraic (b) / Algebraic (2*a);
      ++oi;
    }
    else if (sign_disc > 0)
    {
      // We have two distinct real roots. We return them in ascending order.
      const Algebraic      sqrt_disc = leda::sqrt (Algebraic (disc));
      const Algebraic      alg_b = b;
      const Algebraic      alg_2a = 2*a;

      if (sign_a > 0)
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
    // Compute the true degree for the polynomial.
    while (degree > 0 && leda::sign (coeffs[degree]) == 0)
      degree--;

    LEDA_VECTOR<Integer>  cp (degree + 1);
    unsigned int    i;

    for (i = 0; i <= degree; i++ )
      cp[i] = coeffs[i];
    
    Polynomial            poly (cp);
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
    while (leda::sign (coeffs[degree]) == 0)
    {
      if (degree == 0)
        return (false);
      degree--;
    }
    
    // Compute the least common multiplicand (LCM) of the denominators,
    // denoted L, and the greatest common divisor (GCD) of the numerators,
    // denoted G.
    unsigned int   index = degree;
    Integer        denom_lcm, temp_lcm;
    Integer        numer_gcd, temp_gcd;
    Integer        numer, denom;
    
    denom_lcm = coeffs[index].denominator();
    numer_gcd = coeffs[index].numerator();
    
    index--;
    while (index >= 0)
    {
      numer = coeffs[index].numerator();
      
      if (leda::sign (numer) != 0)
      {
        denom = coeffs[index].denominator();
        
        temp_lcm = denom_lcm;
        temp_gcd = numer_gcd;
        
        denom_lcm *= denom;
        denom_lcm /= leda::gcd (temp_lcm, denom);
        
        numer_gcd = leda::gcd (temp_gcd, numer);
      }
      
      index--;
    }
    
    // Generate the output coefficients (n(i)*L) / (d(i)*G).
    Integer                *z_coeffs = new Integer [degree + 1];
    
    for (index = 0; index <= degree; index++)
    {
      z_coeffs[index] = (coeffs[index].numerator() * denom_lcm) /
                        (numer_gcd * coeffs[index].denominator());
    }
    
    // Set the output.
    poly = construct_polynomial (z_coeffs, degree);
    poly_denom = numer_gcd;
    
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
    // Compute the true degrees for both polynomials.
    unsigned int p_deg = p_degree;
    while (p_deg >= 0 && leda::sign (p_coeffs[p_deg]) == 0)
      p_deg--;

    unsigned int q_deg = q_degree;
    while (q_deg >= 0 && leda::sign (q_coeffs[q_deg]) == 0)
      q_deg--;

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

    // Construct two rational polynomials.
    LEDA_VECTOR<Rational> cp_rat(p_deg + 1);
    LEDA_VECTOR<Rational> cq_rat(q_deg + 1);
    unsigned int                   i;

    for (i = 0; i <= p_deg; i++)
      cp_rat[i] = p_coeffs[i];

    for (i = 0; i <= q_deg; i++)
      cq_rat[i] = q_coeffs[i];

    leda::polynomial<Rational>    P (cp_rat);
    leda::polynomial<Rational>    Q (cq_rat);

    // Compute the GCD of the two polynomials and normalize them.
    leda::polynomial<Rational>    g = leda::poly_gcd (P, Q);
    
    if (g.degree() > 0)
    {
      P = P / g;                   // ???
      p_deg -= g.degree();
      Q = Q / g;                   // ???
      q_deg -= g.degree();
    }

    // Construct two polynomials with integer coefficients.
    Rational      *p_ncfs = new Rational [p_deg + 1];
    Rational      *q_ncfs = new Rational [q_deg + 1];
    Integer        p_scale, q_scale;

    for (i = 0; i <= p_deg; i++)
      p_ncfs[i] = P[i];

    construct_polynomial (p_ncfs, p_deg,
			  p_poly, q_scale);

    for (i = 0; i <= q_deg; i++)
      q_ncfs[i] = Q[i];

    construct_polynomial (q_ncfs, q_deg,
			  q_poly, p_scale);

    delete[] p_ncfs;
    delete[] q_ncfs;

    // Scale the result polynomials.
    p_poly.scale_up (p_scale);
    q_poly.scale_up (q_scale);
   
    return (true);
  }

  /*!
   * Compute the degree of a polynomial.
   */
  int degree (const Polynomial& poly) const
  {
    return (poly.degree());
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

    return (poly[i]);         // RWRW
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
    return (poly.evaluate (x));
  }

  /*!
   * Compute the derivative of the given polynomial.
   * \param poly The polynomial p(x).
   * \return The derivative p'(x).
   */
  Polynomial derive (const Polynomial& poly) const
  {
    return (leda::diff(poly));
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
    poly.scale_up(factor);
    return;
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
      Polynomial q;
      polyA.euclidean_division(polyB, q, rem);
      return q;
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
    int        degree = poly.degree();

    if (degree <= 0)
      return (oi);

    // Check if we really have a simple quadratic equation.
    if (degree <= 2)
    {
      return (solve_quadratic_equation ((degree == 2 ? poly[2] : 0), 
                                        poly[1],
                                        poly[0],
                                        oi));
    }

    // Check if the polynomial is square free. 
    Polynomial             g = leda::poly_gcd(poly, leda::diff(poly));
    const bool             is_square_free = (g.degree() == 0);

    // Compute the real roots of the polynomial.
    leda::list<Algebraic>  roots;
    
    leda::real_roots (poly, roots,
                      leda::real::isolating_algorithm,
                      is_square_free);

    while (! roots.empty())
    {
      *oi = roots.pop_front();
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
    // Check if the polynomial is square free. 
    Polynomial             g = leda::poly_gcd(poly, leda::diff(poly));
    const bool             is_square_free = (g.degree() == 0);
    const Algebraic        x_l = x_min;
    const Algebraic        x_r = x_max;

    // Compute the real roots of the polynomial.
    leda::list<Algebraic>  roots;
    
    leda::real_roots (poly, roots,
                      leda::real::isolating_algorithm,
                      is_square_free);

    while (! roots.empty())
    {
      const Algebraic&      x = roots.pop_front();

      if (x < x_l)
        continue;

      if (x > x_r)
        break;

      *oi = x;
      ++oi;
    }

    return (oi);
  }
};

CGAL_END_NAMESPACE

#endif
