// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/predicates_on_points_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

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
