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
#include <map>
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

  // Type definition for the intersection-point mapping.
  typedef size_t                                  Curve_id;
  typedef std::pair<Curve_id, Curve_id>           Curve_pair;
  typedef std::pair<Point_2, unsigned int>        Intersection_point_2;
  typedef std::list<Intersection_point_2>         Intersection_list;
  typedef std::pair<Intersection_list, bool>      Intersection_info;
  typedef typename Intersection_list::iterator    Intersection_iter;

  /*! \struct
   * An auxiliary functor for comparing pair of curve IDs.
   */
  struct Less_curve_pair
  {
    bool operator() (const Curve_pair& cp1, const Curve_pair& cp2) const
    {
      // Compare the pairs of IDs lexicographically.
      return (cp1.first < cp2.first ||
              (cp1.first == cp2.first && cp1.second < cp2.second));
    }
  };

  typedef std::map<Curve_pair,
                   Intersection_info,
                   Less_curve_pair>               Intersection_map;
  typedef typename Intersection_map::value_type   Intersection_map_entry;
  typedef typename Intersection_map::iterator     Intersection_map_iterator;

private:

  typedef typename Alg_kernel::Point_2            Alg_point_2;
  typedef typename Rat_kernel::Point_2            Rat_point_2;

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
  bool              _inc_to_right; // Is the parameter value increases when
                                   // traversing the subcurve from left to
                                   // right.
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

    // Check if the vaule of the parameter t increase when we traverse the
    // curve from left to right: If the curve is directed to the right, we
    // check if t_src < t_trg, otherwise we check whether t_src > t_trg.
    const Comparison_result    t_res = CGAL::compare (t_src, t_trg);

    CGAL_precondition (t_res != EQUAL);

    if (_dir_right)
      _inc_to_right = (t_res == SMALLER);
    else
      _inc_to_right = (t_res == LARGER);
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
   * \param inter_map Maps curve pairs to lists of their intersection points.
   * \pre p is in the x-range of the arc.
   * \return SMALLER if the point is below the arc;
   *         LARGER if the point is above the arc;
   *         EQUAL if p lies on the arc.
   */
  Comparison_result point_position (const Point_2& p,
                                    Intersection_map& inter_map) const
  {
    // First check whether p has the same x-coordinate as one of the endpoints.
    const Comparison_result  res1 = CGAL::compare (p.x(), _ps.x());

    if (res1 == EQUAL)
      return (CGAL::compare (p.y(), _ps.y()));

    const Comparison_result  res2 = CGAL::compare (p.x(), _pt.x());

    if (res2 == EQUAL)
      return (CGAL::compare (p.y(), _pt.y()));

    // Make sure that p is in the x-range of our subcurve.
    CGAL_precondition (res1 != res2);

    // Check if the supporting curve is one of p's originators. If so, p
    // obviously lies on our subcurve.
    Algebraic        t;

    if (p.is_originator (_curve, t))
    {
      if (_is_in_range (t))
        return (EQUAL);
    }
    else
    {
      // Get of of p's originating curves and compute its intersections
      // with our x-monotone curve.
      CGAL_assertion (p.originators_begin() != p.originators_end());

      typename Point_2::Originator   orig = *(p.originators_begin());
      bool                           do_overlap;
      std::pair<Intersection_iter, Intersection_iter>
        inter_range = _get_curve_intersections (_curve, orig.first,
                                                inter_map, do_overlap);
      
      if (do_overlap)
        return (EQUAL);
      
      // Go over the intersection points and look for p there.
      Intersection_iter    iit;
      Algebraic            s;
      bool                 is_orig;
      
      for (iit = inter_range.first; iit != inter_range.second; ++iit)
      {
        const Point_2&       pt = iit->first;
        
        is_orig = pt.is_originator (orig.first, t);
        CGAL_assertion (is_orig);
        
        if (CGAL::compare (t, orig.second) == EQUAL)
        {
          // Check if this intersection point lies on the parameter range
          // of our subcurve. If so, p lies on the subcurve.
          is_orig = pt.is_originator (_curve, s);
          CGAL_assertion (is_orig);

          if (_is_in_range (s))
            return (EQUAL);
        }
      }
    }

    // If we reached here, p does not lie on the subcurve.
    // TODO: More careful work here ...
    double     app_x = CGAL::to_double (p.x());
    const int  denom = 100000;
    const int  numer = static_cast<int> (app_x * denom);
    Algebraic  y0 = _get_y (Rational (numer, denom));
    
    return (CGAL::compare (p.y(), y0));
  }

  /*!
   * Compare the relative y-position of two x-monotone subcurve to the right
   * of their intersection point.
   * \param cv The other subcurve.
   * \param p The intersection point.
   * \param inter_map Maps curve pairs to lists of their intersection points.
   * \pre p is the common left endpoint of both subcurves.
   * \return SMALLER if (*this) lies below cv to the right of p;
   *         EQUAL in case of an overlap (should not happen);
   *         LARGER if (*this) lies above cv to the right of p.
   */
  Comparison_result compare_to_right (const Self& cv,
                                      const Point_2& p,
                                      Intersection_map& inter_map) const
  {
    CGAL_precondition (CGAL::compare (p.x(), right().x()) != LARGER);
    CGAL_precondition (CGAL::compare (p.x(), cv.right().x()) != LARGER);

    if (this == &cv)
      return (EQUAL);

    // Check for vertical subcurves. A vertical segment is above any other
    // x-monotone subcurve to the right of their common endpoint.
    if (is_vertical())
    {
      if (cv.is_vertical())
        // Both are vertical segments with a common endpoint, so they overlap:
        return (EQUAL);
      
      return (LARGER);
    }
    else if (cv.is_vertical())
    {
      return (SMALLER);
    }

    // Get the parameter value for the point p.
    Nt_traits       nt_traits;
    Algebraic       t0;
    bool            is_orig = p.is_originator (_curve, t0);

    CGAL_assertion (is_orig);

    // Check if both subcurves originate from the same Bezier curve.
    if (_curve.is_same (cv._curve))
    {
      // In this case we know that we have a vertical tangency at t0, so
      // X'(t0) = 0. We evaluate the sign of Y'(t0) in order to find the
      // vertical position of the two subcurves to the right of this point.
      Polynomial   polyY_der = nt_traits.derive (_curve.y_polynomial());
      CGAL::Sign   sign_der = CGAL::sign (nt_traits.evaluate_at (polyY_der,
                                                                 t0));

      CGAL_assertion (sign_der != CGAL::ZERO);

      if (_inc_to_right == cv._inc_to_right)
      {
        std::cout << "cv1 = " << *this << std::endl;
        std::cout << "cv2 = " << cv << std::endl;
        std::cout << "p = (" << p << ") with t0 = " << t0 << std::endl;
      }
      CGAL_assertion (_inc_to_right != cv._inc_to_right);

      if (_inc_to_right)
        return ((sign_der == CGAL::POSITIVE) ? LARGER : SMALLER);
      else
        return ((sign_der == CGAL::NEGATIVE) ? LARGER : SMALLER);
    }
    
    // Get the range of intersection points between the two supporting curves.
    // In particular, check if the curves overlap.
    bool                                            do_overlap = false;
    std::pair<Intersection_iter, Intersection_iter> inter_range;

    inter_range = _get_curve_intersections (_curve, cv._curve,
                                            inter_map,
                                            do_overlap);

    if (do_overlap)
      return (EQUAL);

    // Compare the slopes of the two supporting curves at p. In the general
    // case, the slopes are not equal and their comparison gives us the
    // vertical order to p's right. 
    Comparison_result   slope_res = _compare_slopes (cv, p);

    if (slope_res != EQUAL)
      return (slope_res);

    // Compare the two subcurves by choosing some point to the right of p
    // and compareing the vertical position there.
    Comparison_result   right_res;

    if (CGAL::compare (right().x(), cv.right().x()) != LARGER)
    {
      right_res = _compare_to_side (cv, p,
                                    true,           // Compare to p's right.
                                    inter_range);
      CGAL_assertion (right_res != EQUAL);
    }
    else
    {
      right_res = cv._compare_to_side (*this, p,
                                       true,        // Compare to p's right.
                                       inter_range);
      CGAL_assertion (right_res != EQUAL);
      right_res = ((right_res == LARGER) ? SMALLER : LARGER);
    }

    return (right_res);
  }

  /*!
   * Compare the relative y-position of two x-monotone subcurve to the left
   * of their intersection point.
   * \param cv The other subcurve.
   * \param p The intersection point.
   * \param inter_map Maps curve pairs to lists of their intersection points.
   * \pre p is the common right endpoint of both subcurves.
   * \return SMALLER if (*this) lies below cv to the right of p;
   *         EQUAL in case of an overlap (should not happen);
   *         LARGER if (*this) lies above cv to the right of p.
   */
  Comparison_result compare_to_left (const Self& cv,
                                     const Point_2& p,
                                     Intersection_map& inter_map) const
  {
    CGAL_precondition (CGAL::compare (p.x(), left().x()) != SMALLER);
    CGAL_precondition (CGAL::compare (p.x(), cv.left().x()) != SMALLER);

    if (this == &cv)
      return (EQUAL);
    
    // Check for vertical subcurves. A vertical segment is below any other
    // x-monotone subcurve to the left of their common endpoint.
    if (is_vertical())
    {
      if (cv.is_vertical())
        // Both are vertical segments with a common endpoint, so they overlap:
        return (EQUAL);
      
      return (SMALLER);
    }
    else if (cv.is_vertical())
    {
      return (LARGER);
    }

    // Get the parameter value for the point p.
    Nt_traits       nt_traits;
    Algebraic       t0;
    bool            is_orig = p.is_originator (_curve, t0);

    CGAL_assertion (is_orig);

    // Check if both subcurves originate from the same Bezier curve.
    if (_curve.is_same (cv._curve))
    {
      // In this case we know that we have a vertical tangency at t0, so
      // X'(t0) = 0. We evaluate the sign of Y'(t0) in order to find the
      // vertical position of the two subcurves to the right of this point.
      Polynomial   polyY_der = nt_traits.derive (_curve.y_polynomial());
      CGAL::Sign   sign_der = CGAL::sign (nt_traits.evaluate_at (polyY_der,
                                                                 t0));

      CGAL_assertion (sign_der != CGAL::ZERO);
      CGAL_assertion (_inc_to_right != cv._inc_to_right);

      if (_inc_to_right)
        return ((sign_der == CGAL::NEGATIVE) ? LARGER : SMALLER);
      else
        return ((sign_der == CGAL::POSITIVE) ? LARGER : SMALLER);
    }
    
    // Get the range of intersection points between the two supporting curves.
    // In particular, check if the curves overlap.
    bool                                            do_overlap = false;
    std::pair<Intersection_iter, Intersection_iter> inter_range;

    inter_range = _get_curve_intersections (_curve, cv._curve,
                                            inter_map,
                                            do_overlap);

    if (do_overlap)
      return (EQUAL);

    // Compare the slopes of the two supporting curves at p. In the general
    // case, the slopes are not equal and their comparison gives us the
    // vertical order to p's right; note that we swap the order of the curves
    // to obtains their position to the left.
    Comparison_result   slope_res = cv._compare_slopes (*this, p);

    if (slope_res != EQUAL)
      return (slope_res);

    // Compare the two subcurves by choosing some point to the left of p
    // and compareing the vertical position there.
    Comparison_result   left_res;

    if (CGAL::compare (left().x(), cv.left().x()) != SMALLER)
    {
      left_res = _compare_to_side (cv, p,
                                   false,          // Compare to p's left.
                                   inter_range);
      CGAL_assertion (left_res != EQUAL);
    }
    else
    {
      left_res = cv._compare_to_side (*this, p,
                                      false,       // Compare to p's left.
                                      inter_range);
      CGAL_assertion (left_res != EQUAL);
      left_res = ((left_res == LARGER) ? SMALLER : LARGER);
    }

    return (left_res);
  }

  /*!
   * Check whether the two subcurves are equal (have the same graph).
   * \param cv The other subcurve.
   * \param inter_map Maps curve pairs to lists of their intersection points.
   * \return (true) if the two subcurves have the same graph;
   *         (false) otherwise.
   */
  bool equals (const Self& cv,
               Intersection_map& inter_map) const
  {
    // Check if the two subcurve have overlapping supporting curves.
    if (! _curve.is_same (cv._curve))
    {
      if (! _curve.has_same_support (cv._curve))
        return (false);

      // Mark that the two curves overlap.
      const Curve_id               id1 = _curve.id();
      const Curve_id               id2 = cv._curve.id();
      Curve_pair                   curve_pair;

      if (id1 < id2)
        curve_pair = Curve_pair (id1, id2);
      else
        curve_pair = Curve_pair (id2, id1);
      
      if (inter_map.find (curve_pair) == inter_map.end())
      {
        Intersection_info&  info = inter_map[curve_pair];
        info.second = true;
      }
    }

    // Check for equality of the endpoints.
    return (left().equals (cv.left()) &&
            right().equals (cv.right()));
  }

  /*!
   * Compute the intersections points with the given subcurve.
   * \param cv The other subcurve.
   * \param inter_map Maps curve pairs to lists of their intersection points.
   * \param oi The output iterator.
   * \return The past-the-end iterator.
   */
  template<class OutputIterator>
  OutputIterator intersect (const Self& cv,
                            Intersection_map& inter_map,
                            OutputIterator oi) const
  {
    // In case we have two x-monotone subcurves of the same Bezier curve,
    // the only intersections that may occur are common endpoints.
    if (_curve.is_same (cv._curve))
    {
      if (left().is_same (cv.left()) || left().is_same (cv.right()))
      {
        *oi = CGAL::make_object (Intersection_point_2 (left(), 0));
        ++oi;
      }

      if (right().is_same (cv.left()) || right().is_same (cv.right()))
      {
        *oi = CGAL::make_object (Intersection_point_2 (right(), 0));
        ++oi;
      }

      return (oi);
    }

    // Obtain the intersection points between the supporting Bezier curves.
    bool            do_overlap;
    std::pair<Intersection_iter, Intersection_iter>
      inter_range = _get_curve_intersections (_curve, cv._curve,
                                              inter_map, do_overlap);

    if (do_overlap)
    {
      // TODO: Handle overlapping curves ...
      CGAL_assertion (false);
      return (oi);
    }

    // Report only the intersection points that lie on both given subcurves.
    Intersection_iter    iit;
    Algebraic            s, t;
    bool                 is_orig;

    for (iit = inter_range.first; iit != inter_range.second; ++iit)
    {
      // Get an s-value and a t-value such that _curve(s) == cv._curve(t)
      // is the current intersection point.
      const Point_2&       pt = iit->first;
 
      is_orig = pt.is_originator (_curve, s);
      CGAL_assertion (is_orig);

      is_orig = pt.is_originator (cv._curve, t);
      CGAL_assertion (is_orig);

      // Report on the intersection point only if it is in the range of both
      // subcurves.
      if (_is_in_range (s) && cv._is_in_range (t))
      {
        *oi = CGAL::make_object (*iit);
        ++oi;
      }      
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

    if (! _is_in_range (t0))
      std::cout << *this << std::endl
                << "t0 = " << t0 << std::endl;

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
    // Note that we only allow merging subcurves of the same originating
    // Bezier curve (overlapping curves will not do in this case).
    return (_curve.is_same (cv._curve) &&
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
    CGAL_precondition (_curve.is_same (cv._curve));

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
    // TODO: Is this "legal"? Should we touch the Bezier curve instead
    // so that _trg > _src in all cases?
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
   * Check if the given t-value is in the range of the subcurve.
   */
  bool _is_in_range (const Algebraic& t) const
  {
    const Comparison_result  res1 = CGAL::compare (t, _src);
    const Comparison_result  res2 = CGAL::compare (t, _trg);

    return (res1 == EQUAL || res2 == EQUAL || res1 != res2);
  }

  /*! \struct
   * An auxiliary functor for comparing approximate distances.
   */
  typedef std::pair<double, double>        App_point_2;
  typedef std::pair<Point_2, App_point_2>  Ex_point_2;
  typedef std::list<Ex_point_2>            Point_list;
  typedef typename Point_list::iterator    Point_iter;
  typedef std::pair<double, Point_iter>    Distance_iter;

  struct Less_distance_iter
  {
    bool operator() (const Distance_iter& dit1,
                     const Distance_iter& dit2) const
    {
      return (dit1.first < dit2.first);
    }
  };

  /*!
   * Get the intersection points between two Bezier curves, possibly using
   * a precomputed result.
   * \param B1 The first Bezier curve.
   * \param B2 The second Bezier curve.
   * \param inter_map Maps curve pairs to lists of their intersection points.
   * \param do_overlap Output: Do the two curves overlap.
   * \return A pair of iterators specifying the valid range of intersection
   *         points. 
   */
  std::pair<Intersection_iter, Intersection_iter>
  _get_curve_intersections (const Curve_2& B1,
                            const Curve_2& B2,
                            Intersection_map& inter_map,
                            bool& do_overlap) const
  {
    // Construct the pair of curve IDs. Note that the curve with smaller ID
    // always comes first.
    const Curve_id               id1 = B1.id();
    const Curve_id               id2 = B2.id();
    Curve_pair                   curve_pair;

    if (id1 < id2)
      curve_pair = Curve_pair (id1, id2);
    else
      curve_pair = Curve_pair (id2, id1);
      
    // Try to find the curve pair in the map.
    Intersection_map_iterator    map_iter = inter_map.find (curve_pair);

    if (map_iter == inter_map.end())
    {
      // In case the intersection points between the supporting curves have
      // not been computed before, compute them now and store them in the map.
      Intersection_info&  info = inter_map[curve_pair];
      Intersection_list&  ilist = _intersect_curves (B1, B2, 
                                                     info.first,
                                                     info.second);
      
      do_overlap = info.second;
      return (std::make_pair (ilist.begin(), ilist.end()));
    }

    // Now map_iter points to the desried entry in the map. Use it to obtain
    // the valid range of intersection points.
    do_overlap = map_iter->second.second;
    return (std::make_pair (map_iter->second.first.begin(),
                            map_iter->second.first.end()));
  }

  /*!
   * Compute the intersection points of two Bezier curves.
   * \param B1 The first Bezier curve.
   * \param B2 The second Bezier curve.
   * \param inter_list Output: The list of intersection points (along with
   *                           their multiplicities).
   * \param do_overlap Output: Do the two curves overlap.
   * \return The output list of intersection point (same as inter_list).
   */
  Intersection_list& _intersect_curves (const Curve_2& B1,
                                        const Curve_2& B2,
                                        Intersection_list& inter_list,
                                        bool& do_overlap) const
  {
    inter_list.clear();

    // Compute s-values and t-values such that B1(s) and B2(t) are the
    // intersection points.
    std::list<Algebraic>  s_vals;
    std::list<Algebraic>  t_vals;

    B1.basic_intersect (B2, std::back_inserter (s_vals), do_overlap);

    if (do_overlap)
      // Return an empty list of intersection in case of an overlap.
      return (inter_list);

    B2.basic_intersect (B1, std::back_inserter (t_vals), do_overlap);
    CGAL_assertion (! do_overlap);

    CGAL_assertion (s_vals.size() == t_vals.size());

    // Construct the points according to the s- and t-values. Also compute an
    // approximation for each point.
    typename std::list<Algebraic>::iterator  s_it;
    typename std::list<Algebraic>::iterator  t_it;
    Point_2                                  pt;
    App_point_2                              app_pt;
    Point_list                               pts1;
    Point_list                               pts2;
   
    for (s_it = s_vals.begin(); s_it != s_vals.end(); ++s_it)
    {
      pt = Point_2 (B1, *s_it, false);      // Allow illegal s-values.
      app_pt = App_point_2 (CGAL::to_double (pt.x()),
                            CGAL::to_double (pt.y()));
      pts1.push_back (Ex_point_2 (pt, app_pt));
    }

    for (t_it = t_vals.begin(); t_it != t_vals.end(); ++t_it)
    {
      pt = Point_2 (B2, *t_it, false);      // Allow illegal t-values.
      app_pt = App_point_2 (CGAL::to_double (pt.x()),
                            CGAL::to_double (pt.y()));
      pts2.push_back (Ex_point_2 (pt, app_pt));
    }

    // Go over the points in the pts1 list.
    Point_iter                pit1;
    Point_iter                pit2;
    double                    dx, dy;
    Algebraic                 s, t;
    const Algebraic           one (1);
    bool                      is_orig;
    unsigned int              mult;
    int                       k;

    for (pit1 = pts1.begin(), s_it = s_vals.begin();
         pit1 != pts1.end(); 
         ++pit1, ++s_it)
    {
      // Construct a vector of distances from the current point to all other
      // points in the pts2 list.
      const int                     n_pts2 = pts2.size();
      std::vector<Distance_iter>    dist_vec (n_pts2);

      for (k = 0, pit2 = pts2.begin(); pit2 != pts2.end(); k++, ++pit2)
      {
        // Compute the approximate distance between the teo current points.
        dx = pit1->second.first - pit2->second.first;
        dy = pit1->second.second - pit2->second.second;

        dist_vec[k] = Distance_iter (dx*dx + dy*dy, pit2);
      }
      
      // Sort the vector according to the distances from *pit1.
      std::sort (dist_vec.begin(), dist_vec.end(), 
                 Less_distance_iter());

      // Go over the vector entries, starting from the most distant from *pit1
      // to the closest and eliminate pairs of points (we expect that
      // eliminating the distant points is done easily). We stop when we find
      // a pait for *pit1 or when we are left with a single point.
      Point_2&                p1 = pit1->first;
      bool                    found = false;

      for (k = n_pts2 - 1; !found && k > 0; k--)
      {
        pit2 = dist_vec[k].second;
        const Point_2&          p2 = pit2->first;

        if (p1.equals (p2))
        {
          // Add the originator of *pit2 to the originators list of *pit1, and
          // remove this point from pts2.
          p1.merge_originators (p2);
          pts2.erase (pit2);
          found = true;
        }
      }

      if (! found)
      {
        // We are left with a single point - pair it with *pit1.
        pit2 = dist_vec[0].second;
        p1.merge_originators (pit2->first);
        pts2.erase (pit2);
      }

      // Check if the s- and t-values both lie in the legal range of [0,1].
      // If so, report the updated intersection point.
      is_orig = p1.is_originator (B1, s);
      CGAL_assertion (is_orig);

      if (CGAL::sign (s) != NEGATIVE && CGAL::compare (s, one) != LARGER)
      {
        is_orig = p1.is_originator (B2, t);
        CGAL_assertion (is_orig);
      
        if (CGAL::sign (t) != NEGATIVE && CGAL::compare (t, one) != LARGER)
        {
          // TODO: Currently we give all intersection points multiplicity 0,
          //       stating that we do not know the multiplicity. Can't we at
          //       least identify the points with multiplicity 1?
          mult = 0;
          inter_list.push_back (std::make_pair (p1, mult));
        }
      }
    }

    return (inter_list);
  }

  /*!
   * Compute a y-coordinate of a point on the x-monotone subcurve with a
   * given x-coordinate.
   * \param x0 The given x-coodinate.
   * \return The y-coordinate.
   */
  Algebraic _get_y (const Rational& x0) const
  {
    // Obtain the t-values for with the x-coordinates of the supporting
    // curve equal x0.
    std::list<Algebraic>  t_vals;

    _curve.get_t_at_x (x0, std::back_inserter(t_vals));

    // Find a t-value that is in the range of the current curve.
    Nt_traits                                nt_traits;
    typename std::list<Algebraic>::iterator  t_iter;
    Comparison_result                        res1, res2;

    for (t_iter = t_vals.begin(); t_iter != t_vals.end(); ++t_iter)
    {
      res1 = CGAL::compare (*t_iter, _src);

      if (res1 == EQUAL)
      {
        // Return the y-coordinate of the source point:
        return (_ps.y());
      }

      res2 = CGAL::compare (*t_iter, _trg);

      if (res2 == EQUAL)
      {
        // Return the y-coordinate of the source point:
        return (_pt.y());
      }

      if (res1 != res2)
      {
        // We found a t-value in the range of our x-monotone subcurve.
        // Use this value to compute the y-coordinate.
        return (nt_traits.evaluate_at (_curve.y_polynomial(), *t_iter) /
                nt_traits.convert (_curve.y_norm()));
      }
    }

    // If we reached here, x0 is not in the x-range of our subcurve.
    CGAL_assertion (false);
    return (0);
  }

  /*!
   * Compare the slopes of the subcurve with another given Bezier subcurve at
   * their given intersection point.
   * \param cv The other subcurve.
   * \param p The intersection point.
   * \pre p lies of both subcurves.
   * \pre Neither of the subcurves is a vertical segment.
   * \return SMALLER if (*this) slope is less than cv's;
   *         EQUAL if the two slopes are equal;
   *         LARGER if (*this) slope is greater than cv's.
   */
  Comparison_result _compare_slopes (const Self& cv,
                                     const Point_2& p) const
  {
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
    return (CGAL::compare (slope1, slope2));
  }

  /*!
   * Compare the relative y-position of two x-monotone subcurve to the right
   * (or to the left) of their intersection point, whose multiplicity is
   * greater than 1.
   * \param cv The other subcurve.
   * \param p The query point.
   * \param to_right Should we compare to p's right or to p's left.
   * \param inter_range The range of intersection points between the curves.
   * \pre The x-coordinate of the right endpoint of *this is smaller than
   *      (or equal to) the x-coordinate of the right endpoint of cv.
   * \pre p is the common left endpoint of both subcurves.
   * \pre Neither of the subcurves is a vertical segment.
   * \return SMALLER if (*this) lies below cv to the right of p;
   *         EQUAL in case of an overlap (should not happen);
   *         LARGER if (*this) lies above cv to the right of p.
   */
  Comparison_result _compare_to_side 
    (const Self& cv,
     const Point_2& p,
     bool to_right,
     const std::pair<Intersection_iter, Intersection_iter>& inter_range) const
  {
    // Get the parameter value for the point p.
    Algebraic       t0;
    bool            is_orig = p.is_originator (_curve, t0);

    CGAL_assertion (is_orig);

    // Find the next intersection point that lies to the right of p.
    Intersection_iter    iit;
    Algebraic            next_t, t;
    Comparison_result    res;
    bool                 found = false;

    for (iit = inter_range.first; iit != inter_range.second; ++iit)
    {
      // Check if the current point lies to the right (left) of p. We do so by
      // considering its originating parameter value t.
      const Point_2&       pt = iit->first;

      is_orig = pt.is_originator (_curve, t);
      CGAL_assertion (is_orig);

      res = CGAL::compare (t, t0);
      if ((to_right && ((_inc_to_right && res == LARGER) ||
                        (! _inc_to_right && res == SMALLER))) ||
          (! to_right && ((_inc_to_right && res == SMALLER) ||
                          (! _inc_to_right && res == LARGER))))
      {
        if (! found)
        {
          next_t = t;
          found = true;
        }
        else
        {
          // If we have already located an intersection point to the right
          // (left) of p, choose the leftmost (rightmost) of the two points.
          res = CGAL::compare (t, next_t);
          if ((to_right && ((_inc_to_right && res == SMALLER) ||
                            (! _inc_to_right && res == LARGER))) ||
              (! to_right && ((_inc_to_right && res == LARGER) ||
                              (! _inc_to_right && res == SMALLER))))
          {
            next_t = t;
          }
        }
      }
    }

    // If the next intersection point occurs before the right (left) endpoint
    // of the subcurve, keep it. Otherwise, take the parameter value at
    // the endpoint.
    if (found)
    {
      if (to_right == _dir_right)
        res = CGAL::compare (_trg, next_t);
      else
        res = CGAL::compare (_src, next_t);
    }

    if (! found ||
        (to_right && ((_inc_to_right && res == SMALLER) ||
                      (! _inc_to_right && res == LARGER))) ||
        (! to_right && ((_inc_to_right && res == LARGER) ||
                        (! _inc_to_right && res == SMALLER))))
    {
      next_t = ((to_right == _dir_right) ? _trg : _src);
    }

    // Find a rational value between t0 and t_next. Using this value, we
    // a point with rational coordinates on our subcurve.
    Nt_traits       nt_traits;
    const Rational  mid_t = nt_traits.rational_in_interval (t0, next_t);
    Rat_point_2     q1 = _curve (mid_t);

    // Compute a rational point on *this using the t-value we have computed,
    // and locate a point on cv with the same x-coordinate. We now just have
    // to compare the y-coordinates of the two points.
    Algebraic       y2 = cv._get_y (q1.x());
    
    return (CGAL::compare (nt_traits.convert (q1.y()), y2));
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
