// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_BEZIER_CURVE_2_H
#define CGAL_BEZIER_CURVE_2_H

/*! \file
 * Header file for the _Bezier_curve_2 class.
 */

#include <CGAL/Handle_for.h>
#include <CGAL/Bbox_2.h>
#include <algorithm>
#include <vector>
#include <list>
#include <ostream>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of a Bezier curve, specified by (n+1) control points
 * p_0, ... , p_n that define the curve (X(t), Y(t)) for 0 <= t <= 1,
 * where X(t) and Y(t) are polynomials of degree n.
 *
 * The class is templated with three parameters: 
 * Rat_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of control points (and subsequently also for
 *            the polynomial coefficients). This number type must be the same
 *            as Nt_traits::Rational.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is a number type
 *            for representing algebraic numbers. This number type must be the
 *            same as Nt_traits::Algebraic.
 * Nt_traits A traits class that defines the Integer, Rational and Algebraic
 *           number-types, as well as the Polynomial class (a polynomial with
 *           integer coefficients) and enables performing various operations
 *           on objects of these types.
 */

// Forward declaration:
template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class _Bezier_curve_2;

template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class _Bezier_curve_2_rep
{
  friend class _Bezier_curve_2<Rat_kernel_, Alg_kernel_, Nt_traits_>;

public:

  typedef Rat_kernel_                             Rat_kernel;
  typedef Alg_kernel_                             Alg_kernel;
  typedef Nt_traits_                              Nt_traits;
  
  typedef typename Rat_kernel::Point_2            Rat_point_2;
  typedef typename Alg_kernel::Point_2            Alg_point_2;
  
  typedef typename Nt_traits::Integer             Integer;
  typedef typename Nt_traits::Rational            Rational;
  typedef typename Nt_traits::Algebraic           Algebraic;
  typedef typename Nt_traits::Polynomial          Polynomial;

private:

  // X(t) is given by _polyX(t) / _normX:
  Polynomial        _polyX;       // The polynomial for x.
  Integer           _normX;       // Normalizing factor for y.

  // Y(t) is given by _polyY(t) / _normY:
  Polynomial        _polyY;       // The polynomial for y.
  Integer           _normY;       // Normalizing factor for y.

  typedef std::vector<Rat_point_2>                   Control_point_vec;

  Control_point_vec _ctrl_pts;    // The control points.

  Bbox_2            _bbox;        // A bounding box for the curve.

public:

  /*! Default constructor. */
  _Bezier_curve_2_rep ()
  {}

  /*!
   * Constructor from a given range of control points.
   * \param pts_begin An iterator pointing to the first point in the range.
   * \param pts_end A past-the-end iterator for the range.
   * \pre The value-type of the input iterator must be Rat_kernel::Point_2.
   */
  template <class InputIterator>
  _Bezier_curve_2_rep (InputIterator pts_begin, InputIterator pts_end)
  {
    // Get the number of control points and allocate two vectors for rational
    // coefficients.
    const int        n = std::distance (pts_begin, pts_end) - 1;
    Rational        *coeffsX = new Rational [n + 1];
    Rational        *coeffsY = new Rational [n + 1];
    const Rational   rat_zero = Rational (0);
    int              j, k;

    CGAL_precondition_msg (n > 0,
                           "There must be at least 2 control points.");
 
    _ctrl_pts.resize (n + 1);

    for (j = 0; j <= n; j++)
      coeffsX[j] = coeffsY[j] = rat_zero;

    // Compute the rational coefficients, given by the formula:
    //
    //                     n
    //                   *****
    //                   *   *      / n \   k        n-k
    //   (X(t), Y(t)) =   *    p_k (     ) t  (1 - t)
    //                   *   *      \ k /
    //                   *****
    //                    k=0
    //
    double                 x, y;
    double                 x_min = 0, x_max = 0;
    double                 y_min = 0, y_max = 0;
    Rational               px, py;
    Integer                n_over_k_j;
    bool                   even_exp;

    for (k = 0; pts_begin != pts_end; ++pts_begin, k++)
    {
      // Store the current control point and obtain its coordinates.
      _ctrl_pts[k] = *pts_begin;

      px = pts_begin->x();
      py = pts_begin->y();

      // Update the bounding rectangle of the control points.
      x = CGAL::to_double (px);
      y = CGAL::to_double (py);
      
      if (k == 0)
      {
        x_min = x_max = x;
        y_min = y_max = y;
      }
      else
      {
        if (x < x_min)
          x_min = x;
        else if (x > x_max)
          x_max = x;

        if (y < y_min)
          y_min = y;
        else if (y > y_max)
          y_max = y;
      }

      // By simplifying (1 - t)^(n-k) we obtain that the k'th expression of
      // the above sum is given by:
      //
      //     n-k
      //    *****
      //    *   *     j            n!         j+k
      //     *    (-1)  p_k ---------------- t
      //    *   *            j! k! (n-k-j)!
      //    *****
      //     j=0
      //
      even_exp = true; 
      for (j = 0; j <= n - k; j++)
      {
        n_over_k_j = _choose (n, k, j);

        if (even_exp)
        {
          // We should add the current values to the coefficients of the
          // monomial t^(n_j).
          coeffsX[j + k] += px * n_over_k_j;
          coeffsY[j + k] += py * n_over_k_j;
        }
        else
        {
          // We should subtract the current values from the coefficients of the
          // monomial t^(n_j).
          coeffsX[j + k] -= px * n_over_k_j;
          coeffsY[j + k] -= py * n_over_k_j;
        }

        // As we increment j, negate the "even" flag for the exponent (n-j).
        even_exp = !even_exp;
      } // loop on j.
    } // loop on k.

    // Convert the rational polynomials to polynomials with rational
    // coefficients (plus normalizing factors).
    Nt_traits        nt_traits;

    nt_traits.construct_polynomial (coeffsX, n,
                                    _polyX, _normX);
    delete[] coeffsX;

    nt_traits.construct_polynomial (coeffsY, n,
                                    _polyY, _normY);
    delete[] coeffsY;

    // TODO: Should we normalize the polynomials by GCD computation (?)

    // Update the bounding box.
    _bbox = Bbox_2 (x_min, y_min, x_max, y_max);
  }

private:

  /*!
   * Compute the value of n! / (j! k! (n-k-j)!).
   */
  Integer _choose (int n, int j, int k)
  {
    Integer   reduced_fact = 1;
    Integer   j_fact = 1, k_fact = 1;
    int       i;
    
    for (i = n - k - j + 1; i <= n; i++)
      reduced_fact *= Integer (i);

    for (i = 2; i <= j; i++)
      j_fact *= Integer (i);

    for (i = 2; i <= k; i++)
      k_fact *= Integer (i);

    return (reduced_fact / (j_fact * k_fact));
  }
};

template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class _Bezier_curve_2 :
  public Handle_for<_Bezier_curve_2_rep<Rat_kernel_,
                                        Alg_kernel_,
                                        Nt_traits_> >
{
public:

  typedef Rat_kernel_                             Rat_kernel;
  typedef Alg_kernel_                             Alg_kernel;
  typedef Nt_traits_                              Nt_traits;
  typedef _Bezier_curve_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits>              Self;

private:

  typedef _Bezier_curve_2_rep<Rat_kernel,
                              Alg_kernel,
                              Nt_traits>          Bcv_rep;
  typedef Handle_for<Bcv_rep>                     Bcv_handle;

  typedef typename Bcv_rep::Control_point_vec     Control_pt_vec;

public:

  typedef typename Bcv_rep::Rat_point_2           Rat_point_2;
  typedef typename Bcv_rep::Alg_point_2           Alg_point_2;
  
  typedef typename Bcv_rep::Integer               Integer;
  typedef typename Bcv_rep::Rational              Rational;
  typedef typename Bcv_rep::Algebraic             Algebraic;
  typedef typename Bcv_rep::Polynomial            Polynomial;

  typedef typename Control_pt_vec::const_iterator Control_point_iterator;

public:

  /*!
   * Default constructor.
   */
  _Bezier_curve_2 () :
    Bcv_handle (Bcv_rep())
  {}

  /*!
   * Copy constructor.
   */
  _Bezier_curve_2 (const Self& bc) :
    Bcv_handle (bc)
  {}

  /*!
   * Constructor from a given range of control points.
   * \param pts_begin An iterator pointing to the first point in the range.
   * \param pts_end A past-the-end iterator for the range.
   * \pre The value-type of the input iterator must be Rat_kernel::Point_2.
   */
  template <class InputIterator>
  _Bezier_curve_2 (InputIterator pts_begin, InputIterator pts_end) :
    Bcv_handle (Bcv_rep (pts_begin, pts_end))
  {}

  /*!
   * Assignment operator.
   */
  Self& operator= (const Self& bc)
  {
    if (this == &bc || this->identical (bc))
      return (*this);

    Bcv_handle::operator= (bc);
    return (*this);
  }
  
  /*!
   * Get a unique polynomial ID (based on the actual representation pointer).
   */
  size_t id () const
  {
    const void  *p = reinterpret_cast<const void*> (this->ptr());
    
    return (reinterpret_cast<size_t> (p));
  }

  /*!
   * Get the polynomial for the x-coordinates of the curve.
   */
  const Polynomial& x_polynomial () const
  {
    return (_rep()._polyX);
  }

  /*!
   * Get the normalizing factor for the x-coordinates.
   */
  const Integer& x_norm () const
  {
    return (_rep()._normX);
  }

  /*!
   * Get the polynomial for the y-coordinates of the curve.
   */
  const Polynomial& y_polynomial () const
  {
    return (_rep()._polyY);
  }

  /*!
   * Get the normalizing factor for the y-coordinates.
   */
  const Integer& y_norm () const
  {
    return (_rep()._normY);
  }

  /*!
   * Get the number of control points inducing the Bezier curve.
   */
  unsigned int number_of_control_points () const
  {
    return (_rep()._ctrl_pts.size());
  }

  /*!
   * Get the i'th control point.
   * \pre i must be between 0 and n - 1, where n is the number of control
   *      points.
   */
  const Rat_point_2& control_point (unsigned int i) const
  {
    CGAL_precondition (i < number_of_control_points());

    return ((_rep()._ctrl_pts)[i]);
  }

  /*!
   * Get an interator for the first control point.
   */
  Control_point_iterator control_points_begin () const
  {
    return (_rep()._ctrl_pts.begin());
  }

  /*!
   * Get a past-the-end interator for control points.
   */
  Control_point_iterator control_points_end () const
  {
    return (_rep()._ctrl_pts.end());
  }

  /*!
   * Check if both curve handles refer to the same object.
   */
  bool is_same (const Self& bc) const
  {
    return (this->identical (bc));
  }

  /*!
   * Compute a point of the Bezier curve given a rational t-value.
   * \param t The given t-value.
   * \param check_t Should we check the value of t.
   * \pre If check_t is true, t must be between 0 and 1.
   */
  Rat_point_2 operator() (const Rational& t,
                          bool check_t = true) const
  {
    CGAL_precondition (! check_t ||
                       (CGAL::sign (t) != NEGATIVE &&
                        CGAL::compare (t, Rational(1)) != LARGER));
    
    // Compute the x and y coordinates.
    Nt_traits          nt_traits;
    const Rational     x = nt_traits.evaluate_at (_rep()._polyX, t) /
                           Rational (_rep()._normX, 1);
    const Rational     y = nt_traits.evaluate_at (_rep()._polyY, t) /
                           Rational (_rep()._normY, 1);

    return (Rat_point_2 (x, y));
  }
  
  /*!
   * Compute a point of the Bezier curve given an algebraic t-value.
   * \param t The given t-value.
   * \param check_t Should we check the value of t.
   * \pre If check_t is true, t must be between 0 and 1.
   */
  Alg_point_2 operator() (const Algebraic& t,
                          bool check_t = true) const
  {
    // Check for extermal t values (either 0 or 1).
    Nt_traits          nt_traits;
    const CGAL::Sign   sign_t = CGAL::sign (t);

    CGAL_precondition (! check_t || sign_t != NEGATIVE);

    if (sign_t == ZERO)
    {
      // Is t is 0, simply return the first control point.
      const Rat_point_2&  p_0 = _rep()._ctrl_pts[0];

      return (Alg_point_2 (nt_traits.convert (p_0.x()),
                           nt_traits.convert (p_0.y())));
    }

    Comparison_result  res = CGAL::compare (t, Algebraic(1));

    CGAL_precondition (! check_t || res != LARGER);

    if (res == EQUAL)
    {
      // Is t is 0, simply return the first control point.
      const Rat_point_2&  p_n = _rep()._ctrl_pts[_rep()._ctrl_pts.size() - 1];

      return (Alg_point_2 (nt_traits.convert (p_n.x()),
                           nt_traits.convert (p_n.y())));
    }

    // The t-value is between 0 and 1: Compute the x and y coordinates.
    const Algebraic    x = nt_traits.evaluate_at (_rep()._polyX, t) /
                           nt_traits.convert (_rep()._normX);
    const Algebraic    y = nt_traits.evaluate_at (_rep()._polyY, t) /
                           nt_traits.convert (_rep()._normY);

    return (Alg_point_2 (x, y));
  }
 
  /*!
   * Compute the points with vertical tangents on the curve. The function
   * actually returns t-values such that the tangent at (*this)(t) is vertical.
   * \param oi Output: An output iterator for the t-values.
   * \return A past-the-end iterator for the t-values.
   */
  template <class OutputIterator>
  OutputIterator vertical_tangency_points (OutputIterator oi) const
  {
    // Find all t-values such that X'(t) = 0.
    Nt_traits             nt_traits;
    Polynomial            polyX_der = nt_traits.derive (_rep()._polyX);
    std::list<Algebraic>  t_vals;
    const Algebraic       one = Algebraic(1);

    nt_traits.compute_polynomial_roots (polyX_der, std::back_inserter(t_vals));
    // TODO: Compute polynomial roots in interval [0,1]

    // Take only t-values strictly between 0 and 1. Note that we use the
    // fact that the list of roots we obtain is sorted in ascending order.
    typename std::list<Algebraic>::iterator  t_iter = t_vals.begin();

    while (t_iter != t_vals.end() && CGAL::sign (*t_iter) != POSITIVE)
      ++t_iter;

    while (t_iter != t_vals.end() && CGAL::compare (*t_iter, one) != LARGER)
    {
      *oi = *t_iter;
      ++oi;
      
      ++t_iter;
    }

    return (oi);
  }

  /*!
   * Compute the points with horizontal tangents on the curve. The function
   * actually returns t-values such that the tangent at (*this)(t) is
   * horizontal.
   * \param oi Output: An output iterator for the t-values.
   * \return A past-the-end iterator for the t-values.
   */
  template <class OutputIterator>
  OutputIterator horizontal_tangency_points (OutputIterator oi) const
  {
    // Find all t-values such that Y'(t) = 0.
    Nt_traits             nt_traits;
    Polynomial            polyY_der = nt_traits.derive (_rep()._polyY);
    std::list<Algebraic>  t_vals;
    const Algebraic       one = Algebraic(1);

    nt_traits.compute_polynomial_roots (polyY_der, std::back_inserter(t_vals));
    // TODO: compute roots in [0, 1].

    // Take only t-values strictly between 0 and 1. Note that we use the
    // fact that the list of roots we obtain is sorted in ascending order.
    typename std::list<Algebraic>::iterator  t_iter = t_vals.begin();

    while (t_iter != t_vals.end() && CGAL::sign (*t_iter) != POSITIVE)
      ++t_iter;

    while (t_iter != t_vals.end() && CGAL::compare (*t_iter, one) != LARGER)
    {
      *oi = *t_iter;
      ++oi;
      
      ++t_iter;
    }

    return (oi);
  }

  /*!
   * Compute the intersection points between two Bezier curves. The function
   * returns a list of t-values such that (*this)(t) is an intersection point
   * with the other curve, where t is always in [0, 1].
   * \param bc The other Bezier curve.
   * \param oi Output: An output iterator for the t-values.
   * \param do_ovlp Output: Do the two curves overlap.
   * \return A past-the-end iterator for the t-values.
   */
  template <class OutputIterator>
  OutputIterator intersect (const Self& bc, OutputIterator oi,
                            bool& do_ovlp) const
  {
    // Compute all t-values that correspond to intersection points. 
    std::list<Algebraic>  t_vals;

    basic_intersect (bc, std::back_inserter (t_vals), do_ovlp);

    if (do_ovlp)
      return;

    // Report only t-values in the legal range of [0, 1]. Note that we use
    // the fact that the roots we compute are given in ascending order, so
    // we start from the first non-negative root and stop when we come across
    // a root that is greater than 1.
    const Algebraic                          one = Algebraic(1);
    typename std::list<Algebraic>::iterator  t_iter = t_vals.begin();

    while (t_iter != t_vals.end() && CGAL::sign (*t_iter) == NEGATIVE)
      ++t_iter;

    while (t_iter != t_vals.end() && CGAL::compare (*t_iter, one) != LARGER)
    {
      *oi = *t_iter;
      ++oi;
      
      ++t_iter;
    }

    return (oi);
  }

  /*!
   * Compute the intersection points between two Bezier curves. The function
   * returns a list of t-values such that (*this)(t) is an intersection point
   * with the other curve. Note that this function also returns "imaginary"
   * t-values, namely values outside the range [0, 1].
   * \param bc The other Bezier curve.
   * \param oi Output: An output iterator for the t-values.
   * \param do_ovlp Output: Do the two curves overlap.
   * \return A past-the-end iterator for the t-values.
   */
  template <class OutputIterator>
  OutputIterator basic_intersect (const Self& bc, OutputIterator oi,
                                  bool& do_ovlp) const
  {
    // In case the bounding boxes of the two curves do not overlap, they
    // cannot have any intersection.
    if (! do_overlap (bbox(), bc.bbox()))
    {
      do_ovlp = false;
      return (oi);
    }

    // RWRW: Experimental code:
    if (bbox().xmin() == bc.bbox().xmax() ||
        bbox().xmax() == bc.bbox().xmin() ||
        bbox().ymin() == bc.bbox().ymax() ||
        bbox().ymax() == bc.bbox().ymin())
    {
      // Check if there are equal endpoints.
      Rat_kernel                    ker;
      typename Rat_kernel::Equal_2  equal = ker.equal_2_object();
      const Rat_point_2&  s1 = this->control_point (0);
      const Rat_point_2&  t1 = this->control_point 
                                       (this->number_of_control_points() - 1);
      const Rat_point_2&  s2 = bc.control_point (0);
      const Rat_point_2&  t2 = bc.control_point 
                                       (bc.number_of_control_points() - 1);
      
      if (equal (s1, s2) || equal (s1, t2))
      {
        *oi = Algebraic(0);
        ++oi;
      }
      else if (equal (t1, s2) || equal (t1, t2))
      {
        *oi = Algebraic(1);
        ++oi;
      }
       
      do_ovlp = false;
      return (oi);
    }

    // Let us denote our curve by (X_1(s), Y_1(s)) and the other curve by
    // (X_2(t), Y_2(t)), so we have to solve the system of bivariate
    // polynomials in s and t:
    //    I: X_2(t) - X_1(s) = 0
    //   II: Y_2(t) - Y_1(s) = 0
    //
    Nt_traits                nt_traits;
    Integer                  coeff;
    int                      k;

    // Consruct the bivariate polynomial that corresponds to Equation I:
    const Polynomial&        polyX1 = _rep()._polyX;
    const Integer&           normX1 = _rep()._normX;
    const Polynomial&        polyX2 = bc._rep()._polyX;
    const Integer&           normX2 = bc._rep()._normX;
    const int                degX2 = nt_traits.degree (polyX2);
    std::vector<Polynomial>  coeffsX_st (degX2 + 1);

    for (k = degX2; k >= 0; k--)
    {
      coeff = nt_traits.get_coefficient (polyX2, k) * normX1;
      coeffsX_st[k] = nt_traits.construct_polynomial (&coeff, 0);
    }
    coeffsX_st[0] = coeffsX_st[0] - nt_traits.scale (polyX1, normX2);

    // Consruct the bivariate polynomial that corresponds to Equation II:
    const Polynomial&        polyY1 = _rep()._polyY;
    const Integer&           normY1 = _rep()._normY;
    const Polynomial&        polyY2 = bc._rep()._polyY;
    const Integer&           normY2 = bc._rep()._normY;
    const int                degY2 = nt_traits.degree (polyY2);
    std::vector<Polynomial>  coeffsY_st (degY2 + 1);
    
    for (k = degY2; k >= 0; k--)
    {
      coeff = nt_traits.get_coefficient (polyY2, k) * normY1;
      coeffsY_st[k] = nt_traits.construct_polynomial (&coeff, 0);
    }
    coeffsY_st[0] = coeffsY_st[0] - nt_traits.scale (polyY1, normY2);

    // Compute the resultant of the two bivariate polynomials and obtain
    // a polynomial in s. The report the roots of this polynomial. 
    Polynomial            res = _compute_resultant (coeffsX_st, coeffsY_st);

    if (nt_traits.degree (res) < 0)
    {
      do_ovlp = true;
      return (oi);
    }

    do_ovlp = false;
    return (nt_traits.compute_polynomial_roots (res, oi));
  }

  /*!
   * Print the curve.
   */
  std::ostream& print (std::ostream& os) const
  {    
    // Print the X(t) polynomial.
    os << "(X(t) = ";
    _print_polynomial (os, _rep()._polyX, _rep()._normX, 't');

     // Print the Y(t) polynomial.
    os << ", Y(t) = ";
    _print_polynomial (os, _rep()._polyY, _rep()._normY, 't');
    os << ')';

    return (os);
  }

  /*!
   * Sample a portion of the curve (for drawing purposes, etc.).
   * \param t_start The t-value to start with.
   * \param t_end The t-value to end at.
   * \param n_samples The required number of samples.
   * \param oi Output: An output iterator for the samples. The value-type of
   *                   this iterator must be std::pair<double, double>.
   * \return A past-the-end iterator for the samples.
   */
  template <class OutputIterator>
  OutputIterator sample (const double& t_start, const double& t_end,
                         unsigned int n_samples,
                         OutputIterator oi) const
  {
    // Convert the X(t) and Y(t) polynomial to vectors of double coefficients.
    Nt_traits            nt_traits;
    const int            degX = nt_traits.degree (_rep()._polyX);
    std::vector<double>  coeffsX (degX + 1);
    const double         normX = CGAL::to_double (_rep()._normX);
    int                  k;

    for (k = 0; k <= degX; k++)
      coeffsX[k] = CGAL::to_double (nt_traits.get_coefficient (_rep()._polyX, 
                                                               k)) / normX;

    const int            degY = nt_traits.degree (_rep()._polyY);
    std::vector<double>  coeffsY (degY + 1);
    const double         normY = CGAL::to_double (_rep()._normY);

    for (k = 0; k <= degY; k++)
      coeffsY[k] = CGAL::to_double (nt_traits.get_coefficient (_rep()._polyY, 
                                                               k)) / normY;

    // Sample the approximated curve.
    const int            n = (n_samples >= 2) ? n_samples : 2; 
    const double         delta_t = (t_end - t_start) / (n - 1);
    double               x, y;
    double               t;

    for (k = 0; k < n; k++)
    {
      t = t_start + k * delta_t;
      x = _evaluate_at (coeffsX, degX, t);
      y = _evaluate_at (coeffsY, degY, t);

      *oi = std::make_pair (x, y);
      ++oi;
    }

    return (oi);
  }

  /*!
   * Compute all parameter values t such that the x-coordinate of B(t) is x0.
   * Note that the function does not return only values between 0 and 1, so
   * the output t-values may belong to the imaginary continuation of the curve.
   * \param x0 The given x-coordinate.
   * \param oi Output: An output iterator for the t-values.
   * \return A past-the-end iterator.
   */
  template <class OutputIterator>
  OutputIterator get_t_at_x (const Rational& x0,
                             OutputIterator oi) const
  {
    return (_solve_t_values (_rep()._polyX, _rep()._normX, x0, oi));
  }

  /*!
   * Compute all parameter values t such that the y-coordinate of B(t) is y0.
   * Note that the function does not return only values between 0 and 1, so
   * the output t-values may belong to the imaginary continuation of the curve.
   * \param y0 The given y-coordinate.
   * \param oi Output: An output iterator for the t-values.
   * \return A past-the-end iterator.
   */
  template <class OutputIterator>
  OutputIterator get_t_at_y (const Rational& y0,
                             OutputIterator oi) const
  {
    return (_solve_t_values (_rep()._polyY, _rep._normY(), y0, oi));
  }

  /*!
   * Check if the two curves have the same support.
   */
  bool has_same_support (const Self& bc) const
  {
    // If one curve is of degree d1 and the other of degree d2, there can be
    // at most d1*d2 intersection points between them.
    const int      deg1 = number_of_control_points() - 1;
    const int      deg2 = bc.number_of_control_points() - 1;
    const int      n_samples = deg1*deg2;
    Rat_point_2    p1;
    int            k;

    for (k = 0; k <= n_samples; k++)
    {
      // Compute p1 = B1(k/n_samples), where B1 is (*this) curve.
      if (k == 0)
        p1 = (_rep()._ctrl_pts[0]);
      else if (k == 1)
        p1 = (_rep()._ctrl_pts[_rep()._ctrl_pts.size() - 1]);
      else
        p1 = this->operator() (Rational (k, n_samples));

      // Get all t-values such that the x-coordinate of B2(t) equals x1,
      // and check if there exists a t-value such that the y-coordinate of
      // b2(t) equals the y-coordinate of p1.
      std::list<Algebraic>                           t_vals;
      typename std::list<Algebraic>::const_iterator  t_iter;
      Nt_traits                             nt_traits;
      const Algebraic&                      y1 = nt_traits.convert (p1.y());
      bool                                  eq_y = false;

      bc.get_t_at_x (p1.x(), std::back_inserter(t_vals));

      for (t_iter = t_vals.begin(); t_iter != t_vals.end(); ++t_iter)
      {
        const Alg_point_2&  p2 = bc (*t_iter, false);

        if (CGAL::compare (y1, p2.y()) == CGAL::EQUAL)
        {
          eq_y = true;
          break;
        }
      }
      
      // If we found a point on B1 which is not of B2, the two curves do not
      // have the same support.
      if (! eq_y)
        return (false);
    }
    
    // If we reached here, we found (d1*d2 + 1) common points of B1 and B2.
    // This means they have the same support.
    return (true);
  }

  /*!
   * Get the bounding box of the curve.
   */
  const Bbox_2& bbox () const
  {
    return (_rep()._bbox);
  }

private:

  // Get the representation.
  inline const Bcv_rep& _rep () const
  {
    return (*(this->ptr()));
  }

  inline Bcv_rep& _rep ()
  {
    return (*(this->ptr()));
  }

  /*!
   * Compute the resultant of two bivariate polynomials in u and v with
   * respect to v. The bivariate polynomials are given as vectors of,
   * where bp1[i] is a coefficient of v^i, which is in turn a polynomial in u.
   * \param bp1 The first bivariate polynomial.
   * \param bp2 The second bivariate polynomial.
   * \return The resultant polynomial (a polynomial in u).
   */
  Polynomial _compute_resultant (const std::vector<Polynomial>& bp1,
                                 const std::vector<Polynomial>& bp2) const
  {
    // Create the Sylvester matrix of polynomial coefficients. Also prepare
    // the exp_fact vector, that represents the normalization factor (see
    // below).
    Nt_traits        nt_traits;
    const int        m = bp1.size() - 1;
    const int        n = bp2.size() - 1;
    const int        dim = m + n;
    const Integer    zero = 0;
    const Polynomial zero_poly = nt_traits.construct_polynomial (&zero, 0);
    int              i, j, k;

    std::vector<std::vector<Polynomial> >  mat (dim);
    std::vector <int>                      exp_fact (dim);

    for (i = 0; i < dim; i++)
    {
      mat[i].resize (dim);
      exp_fact[i] = 0;

      for (j = 0; j < dim; j++)
        mat[i][j] = zero_poly;
    }

    // Initialize it with copies of the two bivariate polynomials.
    for (i = 0; i < n; i++)
      for (j = m; j >= 0; j--)
        mat[i][i + j] = bp1[j];

    for (i = 0; i < m; i++)
      for (j = n; j >= 0; j--)
        mat[n + i][i + j] = bp2[j];

    // Perform Gaussian elimination on the Sylvester matrix. The goal is to
    // reach an upper-triangular matrix, whose diagonal elements are mat[0][0]
    // to mat[dim-1][dim-1], such that the determinant of the original matrix
    // is given by:
    //
    //              dim-1
    //             *******
    //              *   *  mat[i][i]
    //              *   *
    //               i=0
    //      ---------------------------------
    //         dim-1
    //        *******            exp_fact[i]
    //         *   *  (mat[i][i])
    //         *   *
    //          i=0
    //
    bool             found_row;
    Polynomial       value;
    int              sign_fact = 1;

    for (i = 0; i < dim; i++)
    {
      // Check if the current diagonal value is a zero polynomial.
      if (nt_traits.degree (mat[i][i]) < 0)
      {
        // If the current diagonal value is a zero polynomial, try to replace
        // the current i'th row with a row with a higher index k, such that
        // mat[k][i] is not a zero polynomial.

        found_row = false;
        for (k = i + 1; k < dim; k++)
        {
          if (nt_traits.degree (mat[k][i]) <= 0)
          {
            found_row = true;
            break;
          }
        }
        
        if (found_row)
        {
          // Swap the i'th and the k'th rows (note that we start from the i'th
          // column, because the first i entries in every row with index i or
          // higher should be zero by now).
          for (j = i; j < dim; j++)
          {
            value = mat[i][j];
            mat[i][j] = mat[k][j];
            mat[k][j] = value;
          }

          // Swapping two rows should change the sign of the determinant.
          // We therefore swap the sign of the normalization factor.
          sign_fact = -sign_fact;
        }
        else
        {
          // In case we could not find a non-zero value, the matrix is
          // singular and its determinant is a zero polynomial.
          return (mat[i][i]);
        }
      }

      // Zero the whole i'th column of the following rows.
      for (k = i + 1; k < dim; k++)
      {
        if (nt_traits.degree (mat[k][i]) >= 0)
        {
          value = mat[k][i];
          mat[k][i] = zero_poly;

          for (j = i + 1; j < dim; j++)
          {
            mat[k][j] = mat[k][j] * mat[i][i] - mat[i][j] * value; 
          }
            
          // We multiplied the current row by the i'th diagonal entry, thus
          // multipling the determinant value by it. We therefore increment
          // the exponent of mat[i][i] in the normalization factor.
          exp_fact[i] = exp_fact[i] + 1;
        }
      }
    }

    // Now, the determinant is simply the product of all diagonal items,
    // divided by the normalizing factor.
    const Integer    sgn (sign_fact);
    Polynomial       det_factor = nt_traits.construct_polynomial (&sgn, 0);
    Polynomial       diag_prod = mat[dim - 1][dim - 1];

    CGAL_assertion (exp_fact [dim - 1] == 0);
    for (i = dim - 2; i >= 0; i--)
    {
      // Try to avoid unnecessary multiplications by ignoring the current
      // diagonal item if its exponent in the normalization factor is greater
      // than 0.
      if (exp_fact[i] > 0)
      {
        exp_fact[i] = exp_fact[i] - 1;
      }
      else
      {
        diag_prod *= mat[i][i];
      }

      for (j = 0; j < exp_fact[i]; j++)
        det_factor *= mat[i][i];
    }

    // In case of a trivial normalization factor, just return the product
    // of diagonal elements.
    if (nt_traits.degree(det_factor) == 0)
      return (diag_prod);

    // Divide the product of diagonal elements by the normalization factor
    // and obtain the determinant (note that we should have no remainder).
    Polynomial       det, rem;

    det = nt_traits.divide (diag_prod, det_factor, rem);
    CGAL_assertion (nt_traits.degree(rem) < 0);
    return (det);
  }

  /*!
   * Compute all parameter values t, such that P(t) = val.
   * \param poly The polynomial.
   * \param norm Its normalizing factor.
   * \param val The required value.
   * \param oi Output: An output iterator for the t-values.
   * \return A past-the-end iterator.
   */
  template <class OutputIterator>
  OutputIterator _solve_t_values (const Polynomial& poly,
                                  const Integer& norm,
                                  const Rational& val,
                                  OutputIterator oi) const
  {
    // Construct the polynomial P(t) - val = 0:
    Nt_traits             nt_traits;
    const Integer         numer = nt_traits.numerator (val);
    const Integer         denom = nt_traits.denominator (val);
    const int             deg = nt_traits.degree (poly);
    Integer              *coeffs = new Integer [deg + 1];
    int                   k;

    for (k = 1; k <= deg; k++)
      coeffs[k] = nt_traits.get_coefficient (poly, k) * denom;
    coeffs[0] = nt_traits.get_coefficient (poly, 0) * denom -
                numer * norm;

    // Solve the polynomial and obtain the t-values.
    OutputIterator  end = nt_traits.compute_polynomial_roots
                            (nt_traits.construct_polynomial (coeffs, deg),
                             oi);

    delete[] coeffs;
    return (end);
  }

  /*!
   * Print a polynomial nicely.
   */
  std::ostream& _print_polynomial (std::ostream& os,
                                   const Polynomial& poly,
                                   const Integer& norm,
                                   char var) const
  {
    Nt_traits   nt_traits;
    const int   deg = nt_traits.degree (poly);
    Rational    coeff;
    CGAL::Sign  sgn;
    int         k;

    if (deg < 0)
    {
      os << '0';
      return (os);
    }

    for (k = deg; k >= 0; k--)
    {
      coeff = Rational (nt_traits.get_coefficient (poly, k), norm);

      if (k == deg)
        os << coeff;
      else if ((sgn = CGAL::sign (coeff)) == POSITIVE)
        os << " + " << coeff;
      else if (sgn == NEGATIVE)
        os << " - " << -coeff;
      else
        continue;

      if (k > 1)
        os << '*' << var << '^' << k;
      else if (k == 1)
        os << '*' << var;
    }

    return (os);
  }

  /*!
   * Evaluate a polynomial with double-precision coefficient at a given point.
   * \param coeffs The coefficients.
   * \param deg The degree of the polynomial.
   * \param t The value to evaluate at.
   * \return The value of the polynomial at t.
   */
  double _evaluate_at (const std::vector<double>& coeffs,
                       int deg, const double& t) const
  {
    // Use Horner's rule to evaluate the polynomial at t.
    double     val = coeffs[deg];
    int        k;

    for (k = deg - 1; k >= 0; k--)
      val = val*t + coeffs[k];

    return (val);
  }
};

/*!
 * Exporter for Bezier curves.
 */
template <class Rat_kernel, class Alg_kernel, class Nt_traits>
std::ostream& 
operator<< (std::ostream& os, 
            const _Bezier_curve_2<Rat_kernel, Alg_kernel, Nt_traits> & bc)
{
  return (bc.print (os));
}

CGAL_END_NAMESPACE

#endif
