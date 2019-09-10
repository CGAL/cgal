// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Ron Wein        <wein@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>
               
#ifndef CGAL_ONE_ROOT_NUMBER_H
#define CGAL_ONE_ROOT_NUMBER_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Header file for the One_root_number<NT> class.
 */
#include <CGAL/Interval_arithmetic.h> 

namespace CGAL {

/*! \class
 * Representation of a one-root number - that is, a number of the form
 * m_alpha + m_beta*sqrt(m_gamma), where m_alpha, m_beta and m_gamma are all rational
 * numbers (actually they should be represented by a number type that supports
 * the operators +, -, * and / in an exact manner).
 */
template <class NumberType_, bool Filter_ = true>
class _One_root_number
{
public:

  typedef NumberType_                       NT;
  typedef _One_root_number<NT, Filter_>     Self;

private:

  // The rational coefficients defining the one-root number:
  NT        m_alpha;
  NT        m_beta;
  NT        m_gamma;
  bool      m_is_rational;      // Is the number rational (that is, m_beta = 0).

public:

  /*! Default constructor. */
  _One_root_number () :
    m_alpha (0),
    m_beta (0),
    m_gamma (0),
    m_is_rational (true)
  {}

  /*! Constructor from a rational number. */
  _One_root_number (const NT& val) :
    m_alpha (val),
    m_beta (0),
    m_gamma (0),
    m_is_rational (true)
  {}

  /*! Construct the one-root number a + b*sqrt(c). */
  _One_root_number (const NT& a, const NT& b, const NT& c) :
    m_alpha (a),
    m_beta (b),
    m_gamma (c),
    m_is_rational (false)
  {
    CGAL_precondition (CGAL::sign(c) == POSITIVE);
  }

  /*! Unary minus. */
  Self operator- () const
  {
    if (m_is_rational)
      return Self (-m_alpha);
    else
      return Self (-m_alpha, -m_beta, m_gamma);
  }

  /*! Add a rational number. */
  Self operator+ (const NT& val) const
  {
    if (m_is_rational)
      return Self (m_alpha + val);
    else
      return Self (m_alpha + val, m_beta, m_gamma);
  }

  /*! Subtract a rational number. */
  Self operator- (const NT& val) const
  {
    if (m_is_rational)
      return Self (m_alpha - val);
    else
      return Self (m_alpha - val, m_beta, m_gamma);
  }

  /*! Multiply by a rational number. */
  Self operator* (const NT& val) const
  {
    if (m_is_rational)
      return Self (m_alpha * val);
    else
      return Self (m_alpha * val, m_beta * val, m_gamma);
  }

  /*! Divide by a rational number. */
  Self operator/ (const NT& val) const
  {
    CGAL_precondition (CGAL::sign(val) != ZERO);

    if (m_is_rational)
      return Self (m_alpha / val);
    else
      return Self (m_alpha / val, m_beta / val, m_gamma);
  }


  /*! Add and assign a rational number. */
  void operator+= (const NT& val)
  {
    m_alpha += val;
    return;
  }

  /*! Subtract and assign a rational number. */
  void operator-= (const NT& val)
  {
    m_alpha -= val;
    return;
  }

  /*! Multiply and assign by a rational number. */
  void operator*= (const NT& val)
  {
    m_alpha *= val;
    if (! m_is_rational)
      m_beta *= val;
    return;
  }

  /*! Divide and assign by a rational number. */
  void operator/= (const NT& val)
  {
    m_alpha /= val;
    if (! m_is_rational)
      m_beta /= val;
    return;
  }

  NT alpha() const
  {
    return (m_alpha);
  }

  NT beta() const
  {
    return (m_beta);
  }

  NT gamma() const
  {
    return (m_gamma);
  }

  bool is_rational() const
  {
    return (m_is_rational);
  }

//private:

  CGAL::Sign _sign () const
  {
    const CGAL::Sign    sign_alpha = CGAL::sign (m_alpha);

    if (m_is_rational)
      return (sign_alpha);

    // If m_alpha and m_beta have the same sign, return this sign.
    const CGAL::Sign    sign_beta = CGAL::sign (m_beta);

    if (sign_alpha == sign_beta)
      return (sign_alpha);

    if (sign_alpha == ZERO)
      return (sign_beta);

    // Compare the squared values of m_alpha and of m_beta*sqrt(m_gamma):
    const Comparison_result  res = CGAL::compare (m_alpha*m_alpha,
                                                  m_beta*m_beta * m_gamma);

    if (res == LARGER)
      return (sign_alpha);
    else if (res == SMALLER)
      return (sign_beta);
    else
      return (ZERO);
  }
};

/*!
 * Compute an isolating interval for the one-root number.
 */
template <class NT, bool FL>
std::pair<double, double> to_interval (const _One_root_number<NT, FL>& x)
{
  const CGAL::Interval_nt<true>   alpha_in = to_interval(x.alpha());
  const CGAL::Interval_nt<true>   beta_in = to_interval(x.beta());
  const CGAL::Interval_nt<true>   gamma_in = to_interval(x.gamma());
  const CGAL::Interval_nt<true>&  x_in = alpha_in + 
                                         (beta_in * CGAL::sqrt(gamma_in));

  return (std::make_pair (x_in.inf(), x_in.sup()));
}

/*!
 * Add a rational number and a one-root number.
 */
template <class NT, bool FL>
_One_root_number<NT, FL> operator+ (const NT& val,
                                    const _One_root_number<NT, FL>& x)
{
  if (x.is_rational())
    return _One_root_number<NT, FL> (val + x.alpha());
  else
    return _One_root_number<NT, FL> (val + x.alpha(), x.beta(), x.gamma());
}

/*!
 * Subtract a one-root number from a rational number.
 */
template <class NT, bool FL>
_One_root_number<NT, FL> operator- (const NT& val,
                                    const _One_root_number<NT, FL>& x)
{
  if (x.is_rational())
    return _One_root_number<NT, FL> (val - x.alpha());
  else
    return _One_root_number<NT, FL> (val - x.alpha(), -x.beta(), x.gamma());
}

/*!
 * Multiply a rational number and a one-root number.
 */
template <class NT, bool FL>
_One_root_number<NT, FL> operator* (const NT& val, 
                                    const _One_root_number<NT, FL>& x)
{
  if (x.is_rational())
    return _One_root_number<NT, FL> (val * x.alpha());
  else
    return _One_root_number<NT, FL> (val * x.alpha(), val * x.beta(), x.gamma());
}

/*!
 * Divide a rational number by a one-root number.
 */
template <class NT, bool FL>
_One_root_number<NT, FL> operator/ (const NT& val,
                                    const _One_root_number<NT, FL>& x)
{
  if (x.is_rational())
  {
    // Simple rational division:
    CGAL_precondition (CGAL::sign (x.alpha()) != ZERO);
    return _One_root_number<NT, FL> (val / x.alpha());
  }

  // Use the fact that:
  //
  //        v            v * (a - b*sqrt(c))
  //   -------------- = ---------------------
  //    a + b*sqrt(c)       a^2 - b^2 * c
  //
  NT   denom = x.alpha()*x.alpha() - x.beta()*x.beta() * x.gamma(); 

  CGAL_precondition (CGAL::sign(denom) != ZERO);

  return _One_root_number<NT, FL> (val * x.alpha() / denom,
                               -val * x.beta() / denom, x.gamma());
}

/*!
 * Get a double-precision approximation of the one-root number.
 */
template <class NT, bool FL>
double to_double (const _One_root_number<NT, FL>& x)
{
  if (x.is_rational())
    return (CGAL::to_double(x.alpha()));

  return (CGAL::to_double(x.alpha()) +
          CGAL::to_double(x.beta()) * std::sqrt(CGAL::to_double(x.gamma()))); 
}

/*!
 * Compute the square of a one-root number.
 */
template <class NT, bool FL>
_One_root_number<NT, FL> square (const _One_root_number<NT, FL>& x)
{
  if (x.is_rational())
    return _One_root_number<NT, FL> (x.alpha() * x.alpha());

  // Use the equation:
  //
  //   (a + b*sqrt(c))^2 = (a^2 + b^2*c) + 2ab*sqrt(c)
  //
  return (_One_root_number<NT, FL> (x.alpha()*x.alpha() + x.beta()*x.beta() * x.gamma(),
                                2 * x.alpha() * x.beta(),
                                x.gamma()));
}

/*!
 * Evaluate the sign of a one-root number.
 */
template <class NT, bool FL>
CGAL::Sign sign (const _One_root_number<NT, FL>& x)
{
  if (FL)
  {
    // Try to filter the sign computation using interval arithmetic.
    const std::pair<double, double>&  x_in = CGAL::to_interval (x); 

    if (x_in.first > 0)
      return (CGAL::POSITIVE);
    else if (x_in.second < 0)
      return (CGAL::NEGATIVE);
  }

  // Perform the exact sign computation.
  return (x._sign());
}

/*!
 * Compare a rational number and a one-root number.
 */
template <class NT, bool FL>
CGAL::Comparison_result compare (const NT& val,
                                 const _One_root_number<NT, FL>& x)
{
  if (x.is_rational())
    return (CGAL::compare (val, x.alpha()));

  if (FL)
  {
    // Try to filter the comparison using interval arithmetic.
    const std::pair<double, double>&  x_in = CGAL::to_interval (val); 
    const std::pair<double, double>&  y_in = CGAL::to_interval (x); 
    
    if (x_in.second < y_in.first)
      return (SMALLER);
    else if (x_in.first > y_in.second)
      return (LARGER);
  }

  // Perform the exact comparison.
  const CGAL::Sign   sgn = (val - x)._sign();

  if (sgn == POSITIVE)
    return (LARGER);
  else if (sgn == NEGATIVE)
    return (SMALLER);
  else
    return (EQUAL);
}

/*!
 * Compare a rational number and a one-root number.
 */
template <class NT, bool FL>
CGAL::Comparison_result compare (const _One_root_number<NT, FL>& x,
                                 const NT& val)
{
  if (x.is_rational())
    return (CGAL::compare (x.alpha(), val));

  if (FL)
  {
    // Try to filter the comparison using interval arithmetic.
    const std::pair<double, double>&  x_in = CGAL::to_interval (x); 
    const std::pair<double, double>&  y_in = CGAL::to_interval (val); 

    if (x_in.second < y_in.first)
      return (SMALLER);
    else if (x_in.first > y_in.second)
      return (LARGER);
  }

  // Perform the exact comparison.
  const CGAL::Sign   sgn = (x - val)._sign();

  if (sgn == POSITIVE)
    return (LARGER);
  else if (sgn == NEGATIVE)
    return (SMALLER);
  else
    return (EQUAL);
}

/*!
 * Compare two one-root numbers.
 */
template <class NT, bool FL>
CGAL::Comparison_result compare (const _One_root_number<NT, FL>& x,
                                 const _One_root_number<NT, FL>& y)
{
  if (x.is_rational())
    return (CGAL::compare (x.alpha(), y));
  else if (y.is_rational())
    return (CGAL::compare (x, y.alpha()));

  if (FL)
  {
    // Try to filter the comparison using interval arithmetic.
    const std::pair<double, double>&  x_in = CGAL::to_interval (x); 
    const std::pair<double, double>&  y_in = CGAL::to_interval (y); 
    
    if (x_in.second < y_in.first)
      return (SMALLER);
    else if (x_in.first > y_in.second)
      return (LARGER);
  }

  // Perform the exact comparison:
  // Note that the comparison of (a1 + b1*sqrt(c1)) and (a2 + b2*sqrt(c2))
  // is equivalent to comparing (a1 - a2) and (b2*sqrt(c2) -  b1*sqrt(c1)).
  // We first determine the signs of these terms.
  const NT          diff_alpha = x.alpha() - y.alpha();
  const CGAL::Sign  sign_left = CGAL::sign (diff_alpha);
  const NT          x_sqr = x.beta()*x.beta() * x.gamma();
  const NT          y_sqr = y.beta()*y.beta() * y.gamma();
  Comparison_result right_res = CGAL::compare (y_sqr, x_sqr);
  CGAL::Sign        sign_right = ZERO;
  
  if (right_res == LARGER)
  {
    // Take the sign of b2:
    sign_right = CGAL::sign (y.beta());
  }
  else if (right_res == SMALLER)
  {
    // Take the opposite sign of b1:
    switch (CGAL::sign (x.beta()))
    {
    case POSITIVE :
      sign_right = NEGATIVE;
      break;
    case NEGATIVE:
      sign_right = POSITIVE;
      break;
    case ZERO:
      sign_right = ZERO;
      break;
    default:
      // We should never reach here.
      CGAL_error();
    }
  }
  else
  {
    // We take the sign of (b2*sqrt(c2) -  b1*sqrt(c1)), where both terms
    // have the same absolute value. The sign is equal to the sign of b2,
    // unless both terms have the same sign, so the whole expression is 0.
    sign_right = CGAL::sign (y.beta());
    if (sign_right == CGAL::sign (x.beta()))
      sign_right = ZERO;
  }

  // Check whether on of the terms is zero. In this case, the comparsion
  // result is simpler:
  if (sign_left == ZERO)
  {
    if (sign_right == POSITIVE)
      return (SMALLER);
    else if (sign_right == NEGATIVE)
      return (LARGER);
    else
      return (EQUAL);
  }
  else if (sign_right == ZERO)
  {
    if (sign_left == POSITIVE)
      return (LARGER);
    else if (sign_left == NEGATIVE)
      return (SMALLER);
    else
      return (EQUAL);
  }

  // If the signs are not equal, we can determine the comparison result:
  if (sign_left != sign_right)
  {
    if (sign_left == POSITIVE)
      return (LARGER);
    else
      return (SMALLER);
  }

  // We now square both terms and look at the sign of the one-root number:
  //   ((a1 - a2)^2 - (b1^2*c1 + b2^2*c2)) + 2*b1*b2*sqrt(c1*c2)
  //
  // If both signs are negative, we should swap the comparsion result
  // we eventually compute.
  const NT          A = diff_alpha*diff_alpha - (x_sqr + y_sqr);
  const NT          B = 2 * x.beta() * y.beta();
  const NT          C = x.gamma() * y.gamma();
  const CGAL::Sign  sgn = (_One_root_number<NT, FL> (A, B, C))._sign();
  const bool        swap_res = (sign_left == NEGATIVE);

  if (sgn == POSITIVE)
    return (swap_res ? SMALLER : LARGER);
  else if (sgn == NEGATIVE)
    return (swap_res ? LARGER : SMALLER);
  else
    return (EQUAL);
}

} //namespace CGAL

#endif
