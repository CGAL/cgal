// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra

// This file's name must begin with a lower-case letter for backward
// compatibility.  Unfortunately, you can't have a file that differs only
// in capitalization on the Windows platforms.

#ifndef CGAL_CONVEX_HULL_CONSTRUCTIVE_TRAITS_2_H
#define CGAL_CONVEX_HULL_CONSTRUCTIVE_TRAITS_2_H

#include <CGAL/license/Convex_hull_2.h>

#include <CGAL/enum.h>

namespace CGAL {
namespace Convex_hulls_2 {
namespace internal {

template <class R>
class Compare_signed_distance_to_cached_line_2
  : public R::Compare_signed_distance_to_line_2
{
  typedef typename R::Compare_signed_distance_to_line_2   Base;

public:
  typedef CGAL::Comparison_result                         result_type;

  typedef typename R::Point_2                             Point;
  typedef typename R::Line_2                              Line;

  Compare_signed_distance_to_cached_line_2(const Base& base_f)
    : Base(base_f),
      line_constructed(false)
  { }

  using Base::operator();

  result_type operator()(const Point& a, const Point& b,
                         const Point& c, const Point& d) const
  {
    if(!line_constructed)
    {
      line_constructed = true;
      l_ab = Line(a,b);
    }

    return operator()(l_ab, c, d);
  }

private:
  mutable bool line_constructed;
  mutable Line l_ab;
};

} // namespace internal
} // namespace Convex_hulls_2

template <class K_>
class Convex_hull_constructive_traits_2
  : public K_
{
public:
  typedef K_                                  K;
  typedef typename K::Point_2                 Point_2;
  typedef typename K::Less_xy_2               Less_xy_2;
  typedef typename K::Less_yx_2               Less_yx_2;
  typedef Convex_hulls_2::internal::Compare_signed_distance_to_cached_line_2<K>
                                              Compare_signed_distance_to_line_2;
  typedef typename K::Less_rotate_ccw_2       Less_rotate_ccw_2;
  typedef typename K::Left_turn_2             Left_turn_2;
  typedef typename K::Equal_2                 Equal_2;
  typedef typename K::Segment_2               Segment_2;

  Less_xy_2
  less_xy_2_object() const
  { return Less_xy_2(); }

  Less_yx_2
  less_yx_2_object() const
  { return Less_yx_2(); }

  Compare_signed_distance_to_line_2
  compare_signed_distance_to_line_2_object() const
  {
    return Compare_signed_distance_to_line_2(
             this->K_::compare_signed_distance_to_line_2_object());
  }

  Less_rotate_ccw_2
  less_rotate_ccw_2_object() const
  { return Less_rotate_ccw_2(); }

  Left_turn_2
  left_turn_2_object() const
  { return Left_turn_2(); }

  Equal_2
  equal_2_object() const
  { return Equal_2(); }
};

// for backward compatibility
template <class K>
class convex_hull_constructive_traits_2
  : public Convex_hull_constructive_traits_2<K>
{ };

} // namespace CGAL

#endif // CGAL_CONVEX_HULL_CONSTRUCTIVE_TRAITS_2_H
