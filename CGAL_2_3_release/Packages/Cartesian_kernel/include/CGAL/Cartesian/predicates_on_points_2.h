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

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
bool
x_equal(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
    return p.x() == q.x();
}

template < class R >
inline
bool
y_equal(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
    return p.y() == q.y();
}

template < class R >
inline
bool
equal_xy(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
    return p.x() == q.x() && p.y() == q.y();
}

template < class R >
inline
bool
less_x(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
    return p.x() < q.x();
}

template < class R >
inline
bool
less_y(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
    return p.y() < q.y();
}

template < class R >
inline
Comparison_result
compare_x(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
    return CGAL_NTS compare(p.x(), q.x());
}

template < class R >
inline
Comparison_result
compare_y(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
    return CGAL_NTS compare(p.y(), q.y());
}

template < class R >
inline
Comparison_result
compare_xy(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_xy(p, q);
}

template < class R >
inline
Comparison_result
compare_deltax_deltay(const PointC2<R CGAL_CTAG>& p,
                      const PointC2<R CGAL_CTAG>& q,
                      const PointC2<R CGAL_CTAG>& r,
                      const PointC2<R CGAL_CTAG>& s)
{
    return compare_deltax_deltayC2(p.x(), q.x(), r.y(), s.y());
}

template < class R >
inline
Comparison_result
compare_lexicographically_xy(const PointC2<R CGAL_CTAG> &p,
                             const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_xyC2(p.x(), p.y(), q.x(), q.y());
}

template < class R >
inline
bool
lexicographically_xy_smaller_or_equal(const PointC2<R CGAL_CTAG> &p,
                                      const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_xy(p, q) != LARGER;
}

template < class R >
inline
bool
lexicographically_xy_larger_or_equal(const PointC2<R CGAL_CTAG> &p,
                                     const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_xy(p, q) != SMALLER;
}

template < class R >
inline
bool
lexicographically_xy_smaller(const PointC2<R CGAL_CTAG> &p,
                             const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_xy(p, q) == SMALLER;
}

template < class R >
inline
bool
lexicographically_xy_larger(const PointC2<R CGAL_CTAG> &p,
                            const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_xy(p, q) == LARGER;
}

template < class R >
inline
Comparison_result
compare_yx(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_xyC2(p.y(), p.x(), q.y(), q.x());
}

template < class R >
inline
Comparison_result
compare_lexicographically_yx(const PointC2<R CGAL_CTAG> &p,
                             const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_xyC2(p.y(), p.x(), q.y(), q.x());
}

template < class R >
inline
bool
lexicographically_yx_smaller_or_equal(const PointC2<R CGAL_CTAG> &p,
                                      const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_yx(p, q) != LARGER;
}

template < class R >
inline
bool
lexicographically_yx_larger_or_equal(const PointC2<R CGAL_CTAG> &p,
                                     const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_yx(p, q) != SMALLER;
}

template < class R >
inline
bool
lexicographically_yx_smaller(const PointC2<R CGAL_CTAG> &p,
                             const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_yx(p, q) == SMALLER;
}

template < class R >
inline
bool
lexicographically_yx_larger(const PointC2<R CGAL_CTAG> &p,
                            const PointC2<R CGAL_CTAG> &q)
{
    return compare_lexicographically_yx(p, q) == LARGER;
}

template < class R >
inline
Orientation
orientation(const PointC2<R CGAL_CTAG> &p,
            const PointC2<R CGAL_CTAG> &q,
            const PointC2<R CGAL_CTAG> &r)
{
    return orientationC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template < class R >
inline
Angle
angle(const PointC2<R CGAL_CTAG> &p,
      const PointC2<R CGAL_CTAG> &q,
      const PointC2<R CGAL_CTAG> &r)
{
    return angleC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template < class R >
inline
bool
collinear(const PointC2<R CGAL_CTAG> &p,
          const PointC2<R CGAL_CTAG> &q,
          const PointC2<R CGAL_CTAG> &r)
{
    return orientation(p, q, r) == COLLINEAR;
}

template < class R >
inline
bool
collinear_are_ordered_along_line(const PointC2<R CGAL_CTAG> &p,
                                 const PointC2<R CGAL_CTAG> &q,
                                 const PointC2<R CGAL_CTAG> &r)
{
    return collinear_are_ordered_along_lineC2
             (p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template < class R >
inline
bool
are_ordered_along_line(const PointC2<R CGAL_CTAG> &p,
                       const PointC2<R CGAL_CTAG> &q,
                       const PointC2<R CGAL_CTAG> &r)
{
    return collinear(p, q, r) && collinear_are_ordered_along_line(p, q, r);
}

template < class R >
inline
bool
collinear_are_strictly_ordered_along_line(const PointC2<R CGAL_CTAG> &p,
                                          const PointC2<R CGAL_CTAG> &q,
                                          const PointC2<R CGAL_CTAG> &r)
{
  return collinear_are_strictly_ordered_along_lineC2
               (p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template < class R >
inline
bool
are_strictly_ordered_along_line(const PointC2<R CGAL_CTAG> &p,
                                const PointC2<R CGAL_CTAG> &q,
                                const PointC2<R CGAL_CTAG> &r)
{
    return collinear(p, q, r) &&
	   collinear_are_strictly_ordered_along_line(p, q, r);
}

template < class R >
inline
bool
left_turn(const PointC2<R CGAL_CTAG> &p,
         const PointC2<R CGAL_CTAG> &q,
         const PointC2<R CGAL_CTAG> &r)
{
    return orientation(p, q, r) == LEFT_TURN;
}

template < class R >
inline
bool
right_turn(const PointC2<R CGAL_CTAG> &p,
          const PointC2<R CGAL_CTAG> &q,
          const PointC2<R CGAL_CTAG> &r)
{
   return orientation(p, q, r) == RIGHT_TURN;
}

template < class R >
inline
bool
left_turn(const Origin &o,
          const PointC2<R CGAL_CTAG> &q,
          const PointC2<R CGAL_CTAG> &r)
{
  return Orientation (sign_of_determinant2x2(qx, qy, rx, ry)) == LEFT_TURN;
}

template < class R >
inline
bool
right_turn(const Origin &o,
           const PointC2<R CGAL_CTAG> &q,
           const PointC2<R CGAL_CTAG> &r)
{
  return Orientation (sign_of_determinant2x2(qx, qy, rx, ry)) == RIGHT_TURN;
}

#ifndef CGAL_NO_DEPRECATED_CODE
template < class R >
inline
bool
leftturn(const PointC2<R CGAL_CTAG> &p,
         const PointC2<R CGAL_CTAG> &q,
         const PointC2<R CGAL_CTAG> &r)
{
    return orientation(p, q, r) == LEFT_TURN;
}

template < class R >
inline
bool
rightturn(const PointC2<R CGAL_CTAG> &p,
          const PointC2<R CGAL_CTAG> &q,
          const PointC2<R CGAL_CTAG> &r)
{
   return orientation(p, q, r) == RIGHT_TURN;
}

template < class R >
inline
bool
leftturn(const Origin &o,
         const PointC2<R CGAL_CTAG> &q,
         const PointC2<R CGAL_CTAG> &r)
{
  return Orientation (sign_of_determinant2x2(qx, qy, rx, ry)) == LEFT_TURN;
}

template < class R >
inline
bool
rightturn(const Origin &o,
          const PointC2<R CGAL_CTAG> &q,
          const PointC2<R CGAL_CTAG> &r)
{
  return Orientation (sign_of_determinant2x2(qx, qy, rx, ry)) == RIGHT_TURN;
}
#endif

template <class R>
inline
Oriented_side
side_of_oriented_circle(const PointC2<R CGAL_CTAG> &p,
                        const PointC2<R CGAL_CTAG> &q,
                        const PointC2<R CGAL_CTAG> &r,
                        const PointC2<R CGAL_CTAG> &t)
{
  return side_of_oriented_circleC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(),
	                           t.x(), t.y());
}

template <class R>
inline
Bounded_side
side_of_bounded_circle(const PointC2<R CGAL_CTAG> &p,
                       const PointC2<R CGAL_CTAG> &q,
                       const PointC2<R CGAL_CTAG> &r,
                       const PointC2<R CGAL_CTAG> &t)
{
  return side_of_bounded_circleC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(),
	                          t.x(), t.y());
}

template <class R>
inline
Bounded_side
side_of_bounded_circle(const PointC2<R CGAL_CTAG> &p,
                       const PointC2<R CGAL_CTAG> &q,
                       const PointC2<R CGAL_CTAG> &t)
{
  return side_of_bounded_circleC2(p.x(), p.y(), q.x(), q.y(), t.x(), t.y());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_POINTS_2_H
