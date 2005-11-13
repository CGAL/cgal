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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_ONE_ROOT_NUMBER_H
#define CGAL_ONE_ROOT_NUMBER_H

/*! \file
 * Header file for the One_root_number<NT> class.
 */

#include <list>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of a one-root number - that is, a number of the form
 * alpha + beta*sqrt(gamma), where alpha, beta and gamma are all rational
 * numbers (actually they should be represented by a number type that supports
 * the operators +, -, * and / in an exact manner).
 */
template <class NumberType_>
class _One_root_number
{
public:

  typedef NumberType_                       NT;
  typedef _One_root_number<NT>              Self;

private:

  // The rational coefficients defining the one-root number:
  NT        alpha;
  NT        beta;
  NT        gamma;
  bool      is_rational;      // Is the number rational (that is, beta = 0).

public:

  /*! Default constructor. */
  _One_root_number () :
    alpha (0),
    beta (0),
    gamma (0),
    is_rational (true)
  {}

  /*! Constructor from a rational number. */
  _One_root_number (const NT& val) :
    alpha (val),
    beta (0),
    gamma (0),
    is_rational (true)
  {}

  /*! Construct the one-root number a + b*sqrt(c). */
  _One_root_number (const NT& a, const NT& b, const NT& c) :
    alpha (a),
    beta (b),
    gamma (c),
    is_rational (false)
  {
    CGAL_precondition (CGAL::sign(c) == POSITIVE);
  }

  /*! Unary minus. */
  Self operator- () const
  {
    if (is_rational)
      return Self (-alpha);
    else
      return Self (-alpha, -beta, gamma);
  }

  /*! Add a rational number. */
  Self operator+ (const NT& val) const
  {
    if (is_rational)
      return Self (alpha + val);
    else
      return Self (alpha + val, beta, gamma);
  }

  /*! Subtract a rational number. */
  Self operator- (const NT& val) const
  {
    if (is_rational)
      return Self (alpha - val);
    else
      return Self (alpha - val, beta, gamma);
  }

  /*! Multiply by a rational number. */
  Self operator* (const NT& val) const
  {
    if (is_rational)
      return Self (alpha * val);
    else
      return Self (alpha * val, beta * val, gamma);
  }

  /*! Divide by a rational number. */
  Self operator/ (const NT& val) const
  {
    CGAL_precondition (CGAL::sign(val) != ZERO);

    if (is_rational)
      return Self (alpha / val);
    else
      return Self (alpha / val, beta / val, gamma);
  }


  /*! Add and assign a rational number. */
  void operator+= (const NT& val)
  {
    alpha += val;
    return;
  }

  /*! Subtract and assign a rational number. */
  void operator-= (const NT& val)
  {
    alpha -= val;
    return;
  }

  /*! Multiply and assign by a rational number. */
  void operator*= (const NT& val)
  {
    alpha *= val;
    if (! is_rational)
      beta *= val;
    return;
  }

  /*! Divide and assign by a rational number. */
  void operator/= (const NT& val)
  {
    alpha /= val;
    if (! is_rational)
      beta /= val;
    return;
  }

  // Friend operators:
  template<class NT_> friend 
  _One_root_number<NT_> operator+ (const NT_& val,
                                   const _One_root_number<NT_>& x);
  
  template<class NT_> friend 
  _One_root_number<NT_> operator- (const NT_& val,
                                   const _One_root_number<NT_>& x);
  
  template<class NT_> friend 
  _One_root_number<NT_> operator* (const NT_& val,
                                   const _One_root_number<NT_>& x);
  
  template<class NT_> friend 
  _One_root_number<NT_> operator/ (const NT_& val,
                                   const _One_root_number<NT_>& x);

  // Friend functions:
  template<class NT_> friend 
  double to_double (const _One_root_number<NT_>& x);

  template<class NT_> friend 
  _One_root_number<NT_> square (const _One_root_number<NT_>& x);

  template<class NT_> friend 
  CGAL::Sign sign (const _One_root_number<NT_>& x);

  template<class NT_> friend 
  CGAL::Comparison_result compare (const NT_& val,
                                   const _One_root_number<NT_>& x);
  
  template<class NT_> friend 
  CGAL::Comparison_result compare (const _One_root_number<NT_>& x,
                                   const NT_& val);

  template<class NT_> friend 
  CGAL::Comparison_result compare (const _One_root_number<NT_>& x,
                                   const _One_root_number<NT_>& y);
};

/*!
 * Add a rational number and a one-root number.
 */
template <class NT>
_One_root_number<NT> operator+ (const NT& val, const _One_root_number<NT>& x)
{
  if (x.is_rational)
    return _One_root_number<NT> (val + x.alpha);
  else
    return _One_root_number<NT> (val + x.alpha, x.beta, x.gamma);
}

/*!
 * Subtract a one-root number from a rational number.
 */
template <class NT>
_One_root_number<NT> operator- (const NT& val, const _One_root_number<NT>& x)
{
  if (x.is_rational)
    return _One_root_number<NT> (val - x.alpha);
  else
    return _One_root_number<NT> (val - x.alpha, -x.beta, x.gamma);
}

/*!
 * Multiply a rational number and a one-root number.
 */
template <class NT>
_One_root_number<NT> operator* (const NT& val, const _One_root_number<NT>& x)
{
  if (x.is_rational)
    return _One_root_number<NT> (val * x.alpha);
  else
    return _One_root_number<NT> (val * x.alpha, val * x.beta, x.gamma);
}

/*!
 * Divide a rational number by a one-root number.
 */
template <class NT>
_One_root_number<NT> operator/ (const NT& val, const _One_root_number<NT>& x)
{
  if (x.is_rational)
  {
    // Simple rational division:
    CGAL_precondition (CGAL::sign (x.alpha) != ZERO);
    return _One_root_number<NT> (val / x.alpha);
  }

  // Use the fact that:
  //
  //        v            v * (a - b*sqrt(c))
  //   -------------- = ---------------------
  //    a + b*sqrt(c)       a^2 - b^2 * c
  //
  NT   denom = x.alpha*x.alpha - x.beta*x.beta * x.gamma; 

  CGAL_precondition (CGAL::sign(denom) != ZERO);

  return _One_root_number<NT> (val * x.alpha / denom,
                               -val * x.beta / denom, x.gamma);
}

/*!
 * Get a double-precision approximation of the one-root number.
 */
template <class NT>
double to_double (const _One_root_number<NT>& x)
{
  if (x.is_rational)
    return (CGAL::to_double(x.alpha));

  return (CGAL::to_double(x.alpha) +
          CGAL::to_double(x.beta) * ::sqrt(CGAL::to_double(x.gamma))); 
}

/*!
 * Compute the square of a one-root number.
 */
template <class NT>
 _One_root_number<NT> square (const _One_root_number<NT>& x)
{
  if (x.is_rational)
    return _One_root_number<NT> (x.alpha * x.alpha);

  // Use the equation:
  //
  //   (a + b*sqrt(c))^2 = (a^2 + b^2*c) + 2ab*sqrt(c)
  //
  const NT          A = x.alpha*x.alpha + x.beta*x.beta * x.gamma;
  const NT          B = 2 * x.alpha * x.beta;

  return (_One_root_number<NT> (x.alpha*x.alpha + x.beta*x.beta * x.gamma,
                                2 * x.alpha * x.beta,
                                x.gamma));
}

/*!
 * Evaluate the sign of a one-root number.
 */
template <class NT>
CGAL::Sign sign (const _One_root_number<NT>& x)
{
  /* RWRW: FILTER
  const double   dx = CGAL::to_double(x);
  if (::fabs (dx) > 0.1)
  {
    if (dx > 0)
      return (POSITIVE);
    else
      return (NEGATIVE);
  }
  */

  const CGAL::Sign    sign_alpha = CGAL::sign (x.alpha);

  if (x.is_rational)
    return (sign_alpha);

  // If alpha and beta have the same sign, return this sign.
  const CGAL::Sign    sign_beta = CGAL::sign (x.beta);

  if (sign_alpha == sign_beta)
    return (sign_alpha);

  if (sign_alpha == ZERO)
    return (sign_beta);

  // Compare the squared values of alpha and of beta*sqrt(gamma):
  const Comparison_result  res = CGAL::compare (x.alpha*x.alpha,
                                                x.beta*x.beta * x.gamma);

  if (res == LARGER)
    return (sign_alpha);
  else if (res == SMALLER)
    return (sign_beta);
  else
    return (ZERO);
}

/*!
 * Compare a rational number and a one-root number.
 */
template <class NT>
CGAL::Comparison_result compare (const NT& val,
                                 const _One_root_number<NT>& x)
{
  if (x.is_rational)
    return (CGAL::compare (val, x.alpha));

  const CGAL::Sign   sgn = CGAL::sign (val - x);

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
template <class NT>
CGAL::Comparison_result compare (const _One_root_number<NT>& x,
                                 const NT& val)
{
  if (x.is_rational)
    return (CGAL::compare (x.alpha, val));

  const CGAL::Sign   sgn = CGAL::sign (x - val);

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
template <class NT>
CGAL::Comparison_result compare (const _One_root_number<NT>& x,
                                 const _One_root_number<NT>& y)
{
  if (x.is_rational)
    return (CGAL::compare (x.alpha, y));
  else if (y.is_rational)
    return (CGAL::compare (x, y.alpha));

  /* RWRW: FILTER
  const double   dx = CGAL::to_double(x);
  const double   dy = CGAL::to_double(y);
  if (::fabs (dx - dy) > 0.1)
  {
    if (dx > dy)
      return (LARGER);
    else
      return (SMALLER);
  }
  */

  // Note that the comparison of (a1 + b1*sqrt(c1)) and (a2 + b2*sqrt(c2))
  // is equivalent to comparing (a1 - a2) and (b2*sqrt(c2) -  b1*sqrt(c1)).
  // We first determine the signs of these terms.
  const NT          diff_alpha = x.alpha - y.alpha;
  const CGAL::Sign  sign_left = CGAL::sign (diff_alpha);
  const NT          x_sqr = x.beta*x.beta * x.gamma;
  const NT          y_sqr = y.beta*y.beta * y.gamma;
  Comparison_result right_res = CGAL::compare (y_sqr, x_sqr);
  CGAL::Sign        sign_right = ZERO;
  
  if (right_res == LARGER)
  {
    // Take the sign of b2:
    sign_right = CGAL::sign (y.beta);
  }
  else if (right_res == SMALLER)
  {
    // Take the opposite sign of b1:
    switch (CGAL::sign (x.beta))
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
      CGAL_assertion (false);
    }
  }
  else
  {
    // We take the sign of (b2*sqrt(c2) -  b1*sqrt(c1)), where both terms
    // have the same absolute value. The sign is equal to the sign of b2,
    // unless both terms have the same sign, so the whole expression is 0.
    sign_right = CGAL::sign (y.beta);
    if (sign_right == CGAL::sign (x.beta))
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
  const NT          B = 2 * x.beta * y.beta;
  const NT          C = x.gamma * y.gamma;
  const CGAL::Sign  sgn = CGAL::sign (_One_root_number<NT> (A, B, C));
  const bool        swap_res = (sign_left == NEGATIVE);

  if (sgn == POSITIVE)
    return (swap_res ? SMALLER : LARGER);
  else if (sgn == NEGATIVE)
    return (swap_res ? LARGER : SMALLER);
  else
    return (EQUAL);
}

CGAL_END_NAMESPACE

#endif
