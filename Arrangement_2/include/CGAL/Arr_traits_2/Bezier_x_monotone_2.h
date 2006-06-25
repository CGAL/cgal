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
 * Header file for the _Bezier_x_monotone_2 class.
 */

#include <CGAL/Arr_trais_2/Bezier_curve_2.h>
#include <CGAL/Arr_trais_2/Bezier_point_2.h>
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

  typedef typename Nt_traits::Algebraic           Algebraic;

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

    // RWRW - To be done!
    return (EQUAL);
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
    // RWRW: To do ...
    return (EQUAL);
  }

  /*!
   * Check whether the two subcurves are equal (have the same graph).
   * \param cv The other subcurve.
   * \return (true) if the two subcurves have the same graph;
   *         (false) otherwise.
   */
  bool equals (const Self& cv) const
  {
    // RWRW: To do ...
    return (false);
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
    // RWRW: To do ...
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
    const bool   is_orig = is_originator (_curve, t0);

    CGAL_precondition (is_orig);
    if (! is_orig)
      return;

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
    // RWRW: to do!
    return (false);
  }

  /*!
   * Merge the current arc with the given arc.
   * \param cv The other subcurve.
   * \pre The two arcs are mergeable.
   */
  void merge (const Self& cv)
  {
    // RWRW: to do!
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
