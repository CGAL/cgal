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

#ifndef CGAL_BEZIER_X_MONOTONE_2_H
#define CGAL_BEZIER_X_MONOTONE_2_H

/*! \file
 * Header file for the _Bezier_x_monotone_2 class.
 */

#include <CGAL/Arr_traits_2/Bezier_curve_2.h>
#include <CGAL/Arr_traits_2/Bezier_point_2.h>
#include <ostream>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Representation of an x-monotone Bezier subcurve, specified by a Bezier curve
 * and t in [t_src, t_trg].
 */
template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class _Bezier_x_monotone_2
{
public:

  typedef Rat_kernel_                             Rat_kernel;
  typedef Alg_kernel_                             Alg_kernel;
  typedef Nt_traits_                              Nt_traits;
  typedef _Bezier_curve_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits>              Curve_2;
  typedef _Bezier_point_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits>              Point_2;
  typedef _Bezier_x_monotone_2<Rat_kernel,
                               Alg_kernel,
                               Nt_traits>         Self;

private:

  typedef typename Nt_traits::Integer             Integer;
  typedef typename Nt_traits::Rational            Rational;
  typedef typename Nt_traits::Algebraic           Algebraic;
  typedef typename Nt_traits::Polynomial          Polynomial;

  // Data members:
  Curve_2           _curve;        // The supporting Bezier curve.
  Algebraic         _src;          // The t-value at the source point.
  Point_2           _ps;           // The source point.
  Algebraic         _trg;          // The t-value at the target point.
  Point_2           _pt;           // The target point.
  bool              _dir_right;    // Is the subcurve directed right (or left).
  bool              _is_vert;      // Is the subcurve a vertical segment.

public:

  /*! Default constructor. */
  _Bezier_x_monotone_2 () :
    _dir_right (false),
    _is_vert (false)
  {}

  /*!
   * Constructor.
   * \param B The supporting Bezier curve.
   * \param t_src The t-value for the source point.
   * \param t_trg The t-value for the target point.
   * \pre The t-values are between 0 and 1.
   */
  _Bezier_x_monotone_2 (const Curve_2& B,
                        const Algebraic& t_src, const Algebraic& t_trg) :
    _curve (B),
    _src (t_src),
    _ps (B, t_src),
    _trg (t_trg),
    _pt (B, t_trg),
    _is_vert (false)
  {
    CGAL_precondition (CGAL::sign (t_src) != NEGATIVE &&
                       CGAL::compare (t_src, Algebraic(1)) != LARGER);
    CGAL_precondition (CGAL::sign (t_trg) != NEGATIVE &&
                       CGAL::compare (t_trg, Algebraic(1)) != LARGER);
    CGAL_precondition (CGAL::compare (t_src, t_trg) != EQUAL);

    // Check if the subcurve is directed left or right.
    const Comparison_result    res = CGAL::compare (_ps.x(), _pt.x());

    if (res == EQUAL)
    {
      // We have a vertical segment. Check if the source is below the target.
      _is_vert = true;
      _dir_right = (CGAL::compare (_ps.y(), _pt.y()) == SMALLER);
    }
    else
    {
      _dir_right = (res == SMALLER);
    }
  }

  /*!
   * Get the supporting Bezier curve.
   */
  const Curve_2& supporting_curve () const
  {
    return (_curve);
  }

  /*!
   * Get the range of t-value over which the subcurve is defined.
   * \return A pair comprised of the t-value for the source point and the
   *         t-value for the target point.
   */
  std::pair<Algebraic, Algebraic> t_range () const
  {
    return (std::make_pair (_src, _trg));
  }

  /*!
   * Get the source point.
   */
  const Point_2& source () const
  {
    return (_ps);
  }

  /*!
   * Get the target point.
   */
  const Point_2& target () const
  {
    return (_pt);
  }

  /*!
   * Get the left endpoint (the lexicographically smaller one).
   */
  const Point_2& left () const
  {
    return (_dir_right ? _ps : _pt);
  }

  /*!
   * Get the right endpoint (the lexicographically larger one).
   */
  const Point_2& right () const
  {
    return (_dir_right ? _pt : _ps);
  }

  /*!
   * Check if the subcurve is a vertical segment.
   */
  bool is_vertical () const
  {
    return (_is_vert);
  }

  /*!
   * Check if the subcurve is directed from left to right.
   */
  bool is_directed_right () const
  {
    return (_dir_right);
  }

  /*!
   * Get the relative position of the query point with respect to the subcurve.
   * \param p The query point.
   * \pre p is in the x-range of the arc.
   * \return SMALLER if the point is below the arc;
   *         LARGER if the point is above the arc;
   *         EQUAL if p lies on the arc.
   */
  Comparison_result point_position (const Point_2& p) const
  {
    CGAL_precondition_code
      (const Comparison_result  res1 = CGAL::compare (p.x(), _ps.x());
       const Comparison_result  res2 = CGAL::compare (p.x(), _pt.x()); );
    CGAL_precondition (res1 == EQUAL || res2 == EQUAL || res1 != res2);

    // RWRW - A hueristic solution. Find a robust one!
    double       t_low, t_high;
    const double init_eps = 0.00000001;
    bool         dir_match;

    if (CGAL::compare (_src, _trg) == SMALLER)
    {
      t_low = CGAL::to_double(_src) - init_eps;
      t_high = CGAL::to_double(_trg) + init_eps;
      dir_match = _dir_right;
    }
    else
    {
      t_low = CGAL::to_double(_trg) - init_eps;
      t_high = CGAL::to_double(_src) + init_eps;
      dir_match = ! _dir_right;
    }

    Nt_traits    nt_traits;
    double       t_mid;
    Algebraic    t_val;
    const double px = CGAL::to_double (p.x());
    Algebraic    x;
    double       diff_x;
    const double eps = 0.0001;
    
    while (true)
    {
      t_mid = (t_low + t_high) / 2;
      t_val = Algebraic (t_mid);
      x = nt_traits.evaluate_at (_curve.x_polynomial(), t_val) /
          nt_traits.convert (_curve.x_norm());

      diff_x = CGAL::to_double(x) - px;
      if (diff_x == 0)
        break;

      if (diff_x > 0)
      {
        if (diff_x < eps)
          break;

        if (dir_match)
          t_high = t_mid;
        else
          t_low = t_mid;
      }
      else
      {
        if (-diff_x < eps)
          break;

       if (dir_match)
          t_low = t_mid;
        else
          t_high = t_mid;
      }
    }

    // Correct the t-value, if necessary.
    if (t_mid < 0)
      t_val = Algebraic(0);
    else if (t_mid > 1)
      t_val = Algebraic(1);

    // Compute the y-coordinate and compare to p.y().
    Algebraic y = nt_traits.evaluate_at (_curve.y_polynomial(), t_val) /
                  nt_traits.convert (_curve.y_norm());

    return (CGAL::compare (p.y(), y));
  }

  /*!
   * Compare the slopes of the subcurve with another given Bezier subcurve at
   * their given intersection point.
   * \param cv The other subcurve.
   * \param p The intersection point.
   * \param mult Output: The mutiplicity of the intersection point.
   * \pre p lies of both subcurves.
   * \pre Neither of the subcurves is a vertical segment.
   * \return SMALLER if (*this) slope is less than cv's;
   *         EQUAL if the two slopes are equal;
   *         LARGER if (*this) slope is greater than cv's.
   */
  Comparison_result compare_slopes (const Self& cv,
				    const Point_2& p,
				    unsigned int& mult) const
  {
    if (_has_same_support (cv))
      return (EQUAL);

    // The slope of (X(t), Y(t)) at t0 is Y'(t)/X'(t).
    // Compute the slope of (*this).
    Algebraic    t1, t2;
    bool         is_orig = p.is_originator (_curve, t1);

    CGAL_assertion (is_orig);

    Nt_traits    nt_traits;
    Polynomial   derivX = nt_traits.derive (_curve.x_polynomial());
    Polynomial   derivY = nt_traits.derive (_curve.y_polynomial());
    Algebraic    slope1 = (nt_traits.evaluate_at (derivY, t1) *
                           nt_traits.convert (_curve.x_norm())) /
                          (nt_traits.evaluate_at (derivX, t1) *
                           nt_traits.convert (_curve.y_norm()));

    // Compute the slope of the other subcurve.
    is_orig = p.is_originator (cv._curve, t2);

    CGAL_assertion (is_orig);

    derivX = nt_traits.derive (cv._curve.x_polynomial());
    derivY = nt_traits.derive (cv._curve.y_polynomial());
    Algebraic    slope2 = (nt_traits.evaluate_at (derivY, t2) *
                           nt_traits.convert (cv._curve.x_norm())) /
                          (nt_traits.evaluate_at (derivX, t2) *
                           nt_traits.convert (cv._curve.y_norm()));

    // Compare the slopes.
    Comparison_result  res = CGAL::compare (slope1, slope2);

    CGAL_assertion (res != EQUAL);  // RWRW: what should we do here?

    mult = 1;
    return (res);
  }

  /*!
   * Check whether the two subcurves are equal (have the same graph).
   * \param cv The other subcurve.
   * \return (true) if the two subcurves have the same graph;
   *         (false) otherwise.
   */
  bool equals (const Self& cv) const
  {
    // Check for equality of the supporting curves and the endpoints.
    return (_has_same_support (cv) &&
            left().equals (cv.left()) &&
            right().equals (cv.right()));
  }

  /*!
   * Compute the intersections points with the given subcurve.
   * \param cv The other subcurve.
   * \param oi The output iterator.
   * \return The past-the-end iterator.
   */
  template<class OutputIterator>
  OutputIterator intersect (const Self& cv,
                            OutputIterator oi) const
  {
    // RWRW: Handle overlapping curves ...

    // Let B_1 be the supporting curve of (*this) and B_2 be the supporting
    // curve of cv. We compute s-values and t-values such that B_1(s) and
    // B_2(t) are the intersection points.
    std::list<Algebraic>  s_vals;
    std::list<Algebraic>  t_vals;

    _curve.intersect (cv._curve, std::back_inserter (s_vals));
    cv._curve.intersect (_curve, std::back_inserter (t_vals));

    CGAL_assertion (s_vals.size() == t_vals.size());

    // Pair the s- and t-values and construct the intersection points.
    typename std::list<Algebraic>::iterator  s_it;
    typename std::list<Algebraic>::iterator  t_it, best_it;
 
    for (s_it = s_vals.begin(); s_it != s_vals.end(); ++s_it)
    {
      if (! _is_in_range (*s_it))
        continue;
      
      // If the point lies in the range of (*this), find a t-value for cv
      // that matches this point.
      Point_2       p (_curve, *s_it);
      const double  px = CGAL::to_double (p.x());
      const double  py = CGAL::to_double (p.y());
      unsigned int  mult = 1;
      double        sqr_dist;
      double        best_sqr_dist = 0;

      best_it = t_vals.end();
      for (t_it = t_vals.begin(); t_it != t_vals.end(); ++t_it)
      {
        Point_2       q (cv._curve, *t_it);
        const double  qx = CGAL::to_double (q.x());
        const double  qy = CGAL::to_double (q.y());

        sqr_dist = (qx - px)*(qx - px) + (qy - py)*(qy - py);
        if (best_it == t_vals.end() || sqr_dist < best_sqr_dist)
        {
          best_it = t_it;
          best_sqr_dist = sqr_dist;
        }
      }

      CGAL_assertion (best_it != t_vals.end());
        
      // If the t-value is in the range of cv, add the intersection point.
      if (cv._is_in_range (*best_it))
      {
        p.add_originator (cv._curve, *best_it);
        *oi = CGAL::make_object (std::make_pair (p, mult));
        ++oi;
      }
      
      // Remove the t-value that matches the current s-value.
      t_vals.erase (best_it);
    }

    return (oi);
  }

  /*!
   * Split the subcurve into two at a given split point.
   * \param p The split point.
   * \param c1 Output: The first resulting arc, lying to the left of p.
   * \param c2 Output: The first resulting arc, lying to the right of p.
   * \pre p lies in the interior of the subcurve (not one of its endpoints).
   */
  void split (const Point_2& p,
              Self& c1, Self& c2) const
  {
    // Duplicate the curve.
    c1 = c2 = *this;
    
    // Find a t-value t0 such that B(t0) is the split point.
    Algebraic    t0;
    const bool   is_orig = p.is_originator (_curve, t0);

    CGAL_precondition (is_orig);
    if (! is_orig)
      return;

    CGAL_precondition (_is_in_range (t0));

    // Perform the split.
    if (_dir_right)
    {
      c1._trg = t0;
      c1._pt = p;
      
      c2._src = t0;
      c2._ps = p;
    }
    else
    {
      c1._src = t0;
      c1._ps = p;
      
      c2._trg = t0;
      c2._pt = p;
    }

    return;
  }

  /*!
   * Check if the two subcurves are mergeable.
   * \param cv The other subcurve.
   * \return Whether the two subcurves can be merged.
   */
  bool can_merge_with (const Self& cv) const
  {
    return (_has_same_support (cv) &&
            (right().equals (cv.left()) || left().equals (cv.right())));
            
    return (false);
  }

  /*!
   * Merge the current arc with the given arc.
   * \param cv The other subcurve.
   * \pre The two arcs are mergeable.
   * \return The merged arc.
   */
  Self merge (const Self& cv) const
  {
    CGAL_precondition (_has_same_support (cv));

    Self    res = cv;

    if (right().equals (cv.left()))
    {
      // Extend the subcurve to the right.
      if (_dir_right)
      {
        res._trg = (cv._dir_right ? cv._trg : cv._src);
        res._pt = cv.right();
      }
      else
      {
        res._src = (cv._dir_right ? cv._trg : cv._src);
        res._ps = cv.right();
      }
    }
    else
    {
      CGAL_precondition (left().equals (cv.right()));

      // Extend the subcurve to the left.
      if (_dir_right)
      {
        res._src = (cv._dir_right ? cv._src : cv._trg);
        res._ps = cv.left();
      }
      else
      {
        res._trg = (cv._dir_right ? cv._src : cv._trg);
        res._pt = cv.left();
      }
    }

    return;
  }

  /*!
   * Flip the subcurve (swap its source and target points).
   * \return The flipped subcurve.
   */
  Self flip () const
  {
    Self  cv = *this;

    cv._src = this->_trg;
    cv._ps = this->_pt;
    cv._trg = this->_src;
    cv._pt = this->_ps;
    cv._dir_right = ! this->_dir_right;

    return (cv);
  }

private:

  /*!
   * Check if the two subcurves have the same supporting curve.
   */
  bool _has_same_support (const Self& cv) const
  {
    return (_curve.equals (cv._curve));
  }

  /*!
   * Check if the given t-value is in the range of the subcurve.
   */
  bool _is_in_range (const Algebraic& t) const
  {
    const Comparison_result  res1 = CGAL::compare (t, _src);
    const Comparison_result  res2 = CGAL::compare (t, _trg);

    return (res1 == EQUAL || res2 == EQUAL || res1 != res2);
  }
};

/*!
 * Exporter for Bezier curves.
 */
template <class Rat_kernel, class Alg_kernel, class Nt_traits>
std::ostream& 
operator<< (std::ostream& os, 
            const _Bezier_x_monotone_2<Rat_kernel, Alg_kernel, Nt_traits>& cv)
{
  std::pair<typename Nt_traits::Algebraic,
            typename Nt_traits::Algebraic>    t_range = cv.t_range();

  os << cv.supporting_curve()
     << " | t in [" << t_range.first 
     << ", " << t_range.second << ']';

  return (os);
}

CGAL_END_NAMESPACE

#endif
