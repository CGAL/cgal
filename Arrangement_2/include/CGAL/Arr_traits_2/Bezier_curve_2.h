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
    Rational               px, py;
    Integer                n_over_k_j;
    bool                   even_exp;

    for (k = 0; pts_begin != pts_end; ++pts_begin, k++)
    {
      px = pts_begin->x();
      py = pts_begin->y();

      // By simplifying (1 - t)^(n-k) we obtain that the k'th expression of
      // the above sum is given by:
      //
      //     n-k
      //    *****
      //    *   *     n-k-j            n!         n-j
      //     *    (-1)      p_k ---------------- t
      //    *   *                j! k! (n-k-j)!
      //    *****
      //     j=0
      //
      even_exp = ((n - k) % 2 == 0); 
      for (j = 0; j <= n - k; j++)
      {
        n_over_k_j = _over (n, k, j);

        if (even_exp)
        {
          // We should add the current values to the coefficients of the
          // monomial t^(n_j).
          coeffsX[n - j] += px * n_over_k_j;
          coeffsY[n - j] += py * n_over_k_j;
        }
        else
        {
          // We should subtract the current values from the coefficients of the
          // monomial t^(n_j).
          coeffsX[n - j] -= px * n_over_k_j;
          coeffsY[n - j] -= py * n_over_k_j;
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
  }

private:

  /*!
   * Compute the value of n! / (j! k! (n-k-j)!).
   */
  Integer _over (int n, int j, int k)
  {
    Integer   reduced_fact = 1;
    Integer   j_fact = 1, k_fact = 1;
    int       i;

    for (i = n - k - j + 1; i <= n; i++)
      reduced_fact *= Integer (i);

    for (i = 1; i < j; i++)
      j_fact += Integer (i);

    for (i = 1; i < k; i++)
      k_fact += Integer (i);

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

public:

  typedef typename Bcv_rep::Rat_point_2           Rat_point_2;
  typedef typename Bcv_rep::Alg_point_2           Alg_point_2;
  
  typedef typename Bcv_rep::Rational              Rational;
  typedef typename Bcv_rep::Algebraic             Algebraic;

private:

  typedef typename Bcv_rep::Integer               Integer;
  typedef typename Bcv_rep::Polynomial            Polynomial;

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
  _Bezier_curve_2 (const Self& bcv) :
    Bcv_handle (bcv)
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
   * Check for equality.
   */
  bool is_same (const Self& bc) const
  {
    return (this->identical (bc));
  }

  /*!
   * Compute a point of the Bezier curve given a t-value.
   * \param t The given t-value.
   * \pre t must be between 0 and 1.
   */
  Alg_point_2 operator() (const Algebraic& t) const
  {
    CGAL_precondition (CGAL::compare (t, Algebraic(0)) != SMALLER &&
                       CGAL::compare (t, Algebraic(1)) != LARGER);

    // Compute the x and y coordinates.
    Nt_traits        nt_traits;
    const Algebraic  x = nt_traits.evaluate_at (_rep()._polyX, t) /
                         nt_traits.convert (_rep()._normX);
    const Algebraic  y = nt_traits.evaluate_at (_rep()._polyY, t) /
                         nt_traits.convert (_rep()._normY);

    return Alg_point_2 (x, y);
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
   * with the other curve.
   * \param bc The other Bezier curve.
   * \param oi Output: An output iterator for the t-values.
   * \return A past-the-end iterator for the t-values.
   */
  template <class OutputIterator>
  OutputIterator intersect (const Self& bc, OutputIterator oi) const
  {
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
    // a polynomial in s. We conider the roots of this resultant polynomial
    // that are between 0 and 1 and obtain the intersection points (note that
    // we use the fact that the roots we compute are given in ascending order).
    Polynomial            res = _compute_resultant (coeffsX_st, coeffsY_st);
    std::list<Algebraic>  s_vals;
    const Algebraic       one = Algebraic(1);

    nt_traits.compute_polynomial_roots (res, std::back_inserter(s_vals));

    typename std::list<Algebraic>::iterator  s_iter = s_vals.begin();

    while (s_iter != s_vals.end() && CGAL::sign (*s_iter) == NEGATIVE)
      ++s_iter;

    while (s_iter != s_vals.end() && CGAL::compare (*s_iter, one) != LARGER)
    {
      *oi = *s_iter;
      ++oi;
      
      ++s_iter;
    }

    return (oi);
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
    // Create the Sylvester matrix of polynomial coefficients.
    Nt_traits        nt_traits;
    const int        m = bp1.size() - 1;
    const int        n = bp2.size() - 1;
    const int        dim = m + n;
    const Integer    zero = 0;
    const Polynomial zero_poly = nt_traits.construct_polynomial (&zero, 0);
    int              i, j, k;

    std::vector<std::vector<Polynomial> >  mat (dim);

    for (i = 0; i < dim; i++)
    {
      mat[i].resize (dim);

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

    // Perform Gaussian elimination on the Sylvester matrix.
    bool             found_row;
    Polynomial       value;
    const Integer    one = 1;
    Polynomial       det_factor = nt_traits.construct_polynomial (&one, 0);
    Polynomial       det_value;

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
          det_factor = -det_factor;
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
          // multipling the determinant value by it.
          det_factor = det_factor * mat[i][i];
        }
      }
    }

    // Now, the determinant is simply the product of all diagonal items,
    // divided by the normalizing factor.
    Polynomial      det, rem;

    det_value = mat[0][0];
    for (i = 1; i < dim; i++)
      det_value = det_value * mat[i][i];

    det = nt_traits.divide (det_value, det_factor, rem);
    CGAL_assertion (nt_traits.degree(rem) < 0);
    return (det);
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
