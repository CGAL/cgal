// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_PREDICATES_ON_POINTS_2_H
#define CGAL_CARTESIAN_PREDICATES_ON_POINTS_2_H

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
bool
x_equal(const PointC2<K> &p, const PointC2<K> &q)
{
  return K().equal_x_2_object()(p, q);
}

template < class K >
inline
bool
y_equal(const PointC2<K> &p, const PointC2<K> &q)
{
  return K().equal_y_2_object()(p, q);
}

template < class K >
inline
bool
equal_xy(const PointC2<K> &p, const PointC2<K> &q)
{
  return p.x() == q.x() && p.y() == q.y();
}

template < class K >
inline
bool
less_x(const PointC2<K> &p, const PointC2<K> &q)
{
  return K().less_x_2_object()(p, q);
}

template < class K >
inline
bool
less_y(const PointC2<K> &p, const PointC2<K> &q)
{
  return K().less_y_2_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_x(const PointC2<K> &p, const PointC2<K> &q)
{
  return K().compare_x_2_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_y(const PointC2<K> &p, const PointC2<K> &q)
{
  return K().compare_y_2_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_xy(const PointC2<K> &p, const PointC2<K> &q)
{
  return K().compare_xy_2_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_deltax_deltay(const PointC2<K>& p,
                      const PointC2<K>& q,
                      const PointC2<K>& r,
                      const PointC2<K>& s)
{
  return compare_deltax_deltayC2(p.x(), q.x(), r.y(), s.y());
}

template < class K >
inline
Comparison_result
compare_lexicographically_xy(const PointC2<K> &p,
                             const PointC2<K> &q)
{
  return K().compare_xy_2_object()(p, q);
}

template < class K >
inline
bool
lexicographically_xy_smaller_or_equal(const PointC2<K> &p,
                                      const PointC2<K> &q)
{
  return compare_lexicographically_xy(p, q) != LARGER;
}

template < class K >
inline
bool
lexicographically_xy_larger_or_equal(const PointC2<K> &p,
                                     const PointC2<K> &q)
{
  return compare_lexicographically_xy(p, q) != SMALLER;
}

template < class K >
inline
bool
lexicographically_xy_smaller(const PointC2<K> &p,
                             const PointC2<K> &q)
{
  return K().less_xy_2_object()(p, q);
}

template < class K >
inline
bool
lexicographically_xy_larger(const PointC2<K> &p,
                            const PointC2<K> &q)
{
  return compare_lexicographically_xy(p, q) == LARGER;
}

template < class K >
inline
Comparison_result
compare_yx(const PointC2<K> &p, const PointC2<K> &q)
{
  return compare_lexicographically_xyC2(p.y(), p.x(), q.y(), q.x());
}

template < class K >
inline
Comparison_result
compare_lexicographically_yx(const PointC2<K> &p,
                             const PointC2<K> &q)
{
  return compare_lexicographically_xyC2(p.y(), p.x(), q.y(), q.x());
}

template < class K >
inline
bool
lexicographically_yx_smaller_or_equal(const PointC2<K> &p,
                                      const PointC2<K> &q)
{
  return compare_lexicographically_yx(p, q) != LARGER;
}

template < class K >
inline
bool
lexicographically_yx_larger_or_equal(const PointC2<K> &p,
                                     const PointC2<K> &q)
{
  return compare_lexicographically_yx(p, q) != SMALLER;
}

template < class K >
inline
bool
lexicographically_yx_smaller(const PointC2<K> &p,
                             const PointC2<K> &q)
{
  return K().less_yx_2_object()(p, q);
}

template < class K >
inline
bool
lexicographically_yx_larger(const PointC2<K> &p,
                            const PointC2<K> &q)
{
  return compare_lexicographically_yx(p, q) == LARGER;
}

template < class K >
inline
Orientation
orientation(const PointC2<K> &p,
            const PointC2<K> &q,
            const PointC2<K> &r)
{
  return K().orientation_2_object()(p, q, r);
}

template < class K >
inline
Angle
angle(const PointC2<K> &p,
      const PointC2<K> &q,
      const PointC2<K> &r)
{
  return K().angle_2_object()(p, q, r);
}

template < class K >
inline
bool
collinear(const PointC2<K> &p,
          const PointC2<K> &q,
          const PointC2<K> &r)
{
  return K().collinear_2_object()(p, q, r);
}

template < class K >
inline
bool
collinear_are_ordered_along_line(const PointC2<K> &p,
                                 const PointC2<K> &q,
                                 const PointC2<K> &r)
{
  return K().collinear_are_ordered_along_line_2_object()(p, q, r);
}

template < class K >
inline
bool
are_ordered_along_line(const PointC2<K> &p,
                       const PointC2<K> &q,
                       const PointC2<K> &r)
{
  return K().are_ordered_along_line_2_object()(p, q, r);
}

template < class K >
inline
bool
collinear_are_strictly_ordered_along_line(const PointC2<K> &p,
                                          const PointC2<K> &q,
                                          const PointC2<K> &r)
{
  return K().collinear_are_strictly_ordered_along_line_2_object()(p, q, r);
}

template < class K >
inline
bool
are_strictly_ordered_along_line(const PointC2<K> &p,
                                const PointC2<K> &q,
                                const PointC2<K> &r)
{
  return K().are_strictly_ordered_along_line_2_object()(p, q, r);
}

template < class K >
inline
bool
left_turn(const PointC2<K> &p,
	  const PointC2<K> &q,
	  const PointC2<K> &r)
{
  return K().left_turn_2_object()(p, q, r);
}

template < class K >
inline
bool
right_turn(const PointC2<K> &p,
	   const PointC2<K> &q,
	   const PointC2<K> &r)
{
  return orientation(p, q, r) == RIGHT_TURN;
}

#ifndef CGAL_NO_DEPRECATED_CODE
template < class K >
inline
bool
leftturn(const PointC2<K> &p,
         const PointC2<K> &q,
         const PointC2<K> &r)
{
  bool THIS_FUNCTION_IS_DEPRECATED; // Use left_turn instead.
  return orientation(p, q, r) == LEFT_TURN;
}

template < class K >
inline
bool
rightturn(const PointC2<K> &p,
          const PointC2<K> &q,
          const PointC2<K> &r)
{
  bool THIS_FUNCTION_IS_DEPRECATED; // Use right_turn instead.
  return orientation(p, q, r) == RIGHT_TURN;
}
#endif

template <class K>
inline
Oriented_side
side_of_oriented_circle(const PointC2<K> &p,
                        const PointC2<K> &q,
                        const PointC2<K> &r,
                        const PointC2<K> &t)
{
  return K().side_of_oriented_circle_2_object()(p, q, r, t);
}

template <class K>
inline
Bounded_side
side_of_bounded_circle(const PointC2<K> &p,
                       const PointC2<K> &q,
                       const PointC2<K> &r,
                       const PointC2<K> &t)
{
  return K().side_of_bounded_circle_2_object()(p, q, r, t);
}

template <class K>
inline
Bounded_side
side_of_bounded_circle(const PointC2<K> &p,
                       const PointC2<K> &q,
                       const PointC2<K> &t)
{
  return K().side_of_bounded_circle_2_object()(p, q, t);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_POINTS_2_H
